'use client';

import React, { useEffect, useRef, useState, useCallback } from 'react';

// ═══════════════════════════════════════════════════════════════════════════
//  CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════

const HELIX_COLOR   = '#E65265';
const HELIX_LIGHT   = '#F0949E';
const HELIX_DARK    = '#B83A4C';
const SHEET_COLOR   = '#F8C630';
const SHEET_LIGHT   = '#FADB6E';
const SHEET_DARK    = '#C49A18';
const LOOP_COLOR    = '#8CB4D1';
const LOOP_LIGHT    = '#B0D0E8';
const LOOP_DARK     = '#5A8FAF';
const HIGHLIGHT_CLR = '#EF4444';
const LIGAND_CLR    = '#F59E0B';
const BG            = '#ffffff';

const CHAIN_PALETTE = [
  '#10B981','#F59E0B','#EF4444','#8B5CF6','#EC4899',
  '#06B6D4','#F97316','#84CC16','#6366F1','#14B8A6',
  '#E11D48','#D97706','#7C3AED','#0891B2','#65A30D',
];

const chainColor = (id: string) => CHAIN_PALETTE[id.charCodeAt(0) % CHAIN_PALETTE.length];

type DisplayMode = 'cartoon' | 'backbone' | 'ballstick';

const MODES: { id: DisplayMode; label: string; icon: string }[] = [
  { id: 'cartoon',   label: 'Cartoon',    icon: '🎞' },
  { id: 'ballstick', label: 'Ball&Stick',  icon: '⚛' },
  { id: 'backbone',  label: 'Backbone',   icon: '🦴' },
];

// ═══════════════════════════════════════════════════════════════════════════
//  TYPES
// ═══════════════════════════════════════════════════════════════════════════

interface HighlightedRes { chain: string; res_num: number; res_name?: string }

interface Props {
  pdbText: string;
  highlightedResidues?: HighlightedRes[];
  selectedChain?: string;
  height?: string;
  className?: string;
  onResidueClick?: (r: { chain: string; res_num: number; res_name: string }) => void;
}

interface Atom {
  x: number; y: number; z: number;
  serial: number; name: string; resName: string;
  chain: string; resSeq: number;
  elem: string; het: boolean; isCA: boolean;
}

type SS = 'helix' | 'sheet' | 'loop';

interface Model {
  atoms: Atom[];
  chains: string[];
  caByChain: Map<string, Atom[]>;
  resByChain: Map<string, Atom[]>;
  ssByChain: Map<string, Map<number, SS>>;
  center: [number, number, number];
  size: number;
}

interface Rot { ax: number; ay: number }

// ═══════════════════════════════════════════════════════════════════════════
//  PDB PARSER
// ═══════════════════════════════════════════════════════════════════════════

function parsePdb(text: string): Model | null {
  const atoms: Atom[] = [];
  const chainSet = new Set<string>();
  const resByChain = new Map<string, Atom[]>();
  const caByChain = new Map<string, Atom[]>();
  const bbNames = new Set(['N','CA','C','O']);
  const helices: { ch: string; s: number; e: number }[] = [];
  const sheets: { ch: string; s: number; e: number }[] = [];

  for (const raw of text.split('\n')) {
    const rec = raw.substring(0, 6).trim();

    if (rec === 'HELIX') {
      const ch = raw.substring(19, 20).trim();
      const s  = parseInt(raw.substring(21, 25));
      const e  = parseInt(raw.substring(33, 37));
      if (ch && !isNaN(s) && !isNaN(e) && e >= s) helices.push({ ch, s, e });
      continue;
    }
    if (rec === 'SHEET') {
      const ch = raw.substring(19, 20).trim();
      const s  = parseInt(raw.substring(21, 25));
      const e  = parseInt(raw.substring(33, 37));
      if (ch && !isNaN(s) && !isNaN(e) && e >= s) sheets.push({ ch, s, e });
      continue;
    }
    if (rec !== 'ATOM' && rec !== 'HETATM') continue;

    const name   = raw.substring(12, 16).trim();
    const resName = raw.substring(17, 20).trim();
    const ch     = raw.substring(21, 22).trim();
    const resSeq = parseInt(raw.substring(22, 26));
    const x      = parseFloat(raw.substring(30, 38));
    const y      = parseFloat(raw.substring(38, 46));
    const z      = parseFloat(raw.substring(46, 54));

    let elem = raw.length >= 78 ? raw.substring(76, 78).trim() : '';
    if (!elem) { for (const c of name) { if (/[A-Za-z]/.test(c)) { elem = c.toUpperCase(); break; } } }

    if (isNaN(x) || isNaN(y) || isNaN(z) || !ch) continue;

    const isCA = name === 'CA';
    const het  = rec === 'HETATM';
    const atom: Atom = { x, y, z, serial: atoms.length, name, resName, chain: ch, resSeq, elem, het, isCA };
    atoms.push(atom);
    chainSet.add(ch);

    if (!het) {
      if (!resByChain.has(ch)) resByChain.set(ch, []);
      resByChain.get(ch)!.push(atom);
      if (isCA) {
        if (!caByChain.has(ch)) caByChain.set(ch, []);
        caByChain.get(ch)!.push(atom);
      }
    }
  }

  if (atoms.length === 0) return null;

  let mnx=Infinity,mxx=-Infinity,mny=Infinity,mxy=-Infinity,mnz=Infinity,mxz=-Infinity;
  for (const a of atoms) {
    if(a.x<mnx)mnx=a.x;if(a.x>mxx)mxx=a.x;
    if(a.y<mny)mny=a.y;if(a.y>mxy)mxy=a.y;
    if(a.z<mnz)mnz=a.z;if(a.z>mxz)mxz=a.z;
  }

  const ssByChain = new Map<string, Map<number, SS>>();
  for (const ch of chainSet) {
    ssByChain.set(ch, assignSS(caByChain.get(ch) || [], helices, sheets));
  }

  return {
    atoms, chains: Array.from(chainSet), caByChain, resByChain, ssByChain,
    center: [(mnx+mxx)/2,(mny+mxy)/2,(mnz+mxz)/2],
    size: Math.max(mxx-mnx, mxy-mny, mxz-mnz, 1),
  };
}

// ═══════════════════════════════════════════════════════════════════════════
//  SECONDARY STRUCTURE
// ═══════════════════════════════════════════════════════════════════════════

function dist(a: Atom, b: Atom) {
  return Math.sqrt((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2);
}

function assignSS(ca: Atom[], helices: { ch:string; s:number; e:number }[], sheets: { ch:string; s:number; e:number }[]): Map<number,SS> {
  const m = new Map<number,SS>();
  for (const a of ca) m.set(a.resSeq, 'loop');

  // Use PDB records if available
  if (helices.length || sheets.length) {
    for (const h of helices) for (const a of ca) if (a.chain===h.ch && a.resSeq>=h.s && a.resSeq<=h.e) m.set(a.resSeq,'helix');
    for (const s of sheets) for (const a of ca) if (a.chain===s.ch && a.resSeq>=s.s && a.resSeq<=s.e) m.set(a.resSeq,'sheet');
    return m;
  }

  // Heuristic detection
  const n = ca.length;
  if (n < 5) return m;

  // Helices: CA(i) to CA(i+4) distance < 7.2 Å
  const flags = new Map<number,number>();
  for (let i = 0; i < n-1; i++) {
    if (ca[i+1].resSeq - ca[i].resSeq !== 1) continue;
    if (i+4 < n && ca[i+4].resSeq - ca[i].resSeq === 4) {
      if (dist(ca[i], ca[i+4]) < 7.2) {
        for (let j=i;j<=i+4&&j<n;j++) flags.set(ca[j].resSeq, (flags.get(ca[j].resSeq)||0)+1);
      }
    }
  }
  for (const [r,c] of flags) if (c>=2) m.set(r,'helix');

  // Clean helix runs (min 4 residues)
  const helSeqs = Array.from(m.entries()).filter(([,t])=>t==='helix').map(([r])=>r).sort((a,b)=>a-b);
  const helRuns: [number,number][] = [];
  if (helSeqs.length) {
    let s=helSeqs[0], e=helSeqs[0];
    for (let i=1;i<helSeqs.length;i++) { if(helSeqs[i]-e<=2)e=helSeqs[i]; else{helRuns.push([s,e]);s=helSeqs[i];e=helSeqs[i];} }
    helRuns.push([s,e]);
  }
  for (const k of m.keys()) if(m.get(k)==='helix') m.set(k,'loop');
  for (const [s,e] of helRuns) if(e-s+1>=4) for(const a of ca) if(a.resSeq>=s&&a.resSeq<=e) m.set(a.resSeq,'helix');

  // Sheets: CA(i) to CA(i+2) distance > 9.5 Å (extended strand)
  for (let i=0;i<n-2;i++) {
    if(ca[i+2].resSeq-ca[i].resSeq!==2) continue;
    if(dist(ca[i],ca[i+2])>9.5 && m.get(ca[i].resSeq)==='loop' && m.get(ca[i+1].resSeq)==='loop' && m.get(ca[i+2].resSeq)==='loop') {
      m.set(ca[i].resSeq,'sheet'); m.set(ca[i+1].resSeq,'sheet'); m.set(ca[i+2].resSeq,'sheet');
    }
  }

  // Clean sheet runs (min 3 residues)
  const shSeqs = Array.from(m.entries()).filter(([,t])=>t==='sheet').map(([r])=>r).sort((a,b)=>a-b);
  const shRuns: [number,number][] = [];
  if (shSeqs.length) {
    let s=shSeqs[0],e=shSeqs[0];
    for(let i=1;i<shSeqs.length;i++){if(shSeqs[i]-e<=2)e=shSeqs[i];else{shRuns.push([s,e]);s=shSeqs[i];e=shSeqs[i];}}
    shRuns.push([s,e]);
  }
  for(const k of m.keys()) if(m.get(k)==='sheet') m.set(k,'loop');
  for(const [s,e] of shRuns) if(e-s+1>=3) for(const a of ca) if(a.resSeq>=s&&a.resSeq<=e) m.set(a.resSeq,'sheet');

  return m;
}

// ═══════════════════════════════════════════════════════════════════════════
//  3D MATH
// ═══════════════════════════════════════════════════════════════════════════

function proj(a: Atom, cx:number,cy:number,cz:number, scale:number, rot:Rot, cw:number, ch:number) {
  let dx=a.x-cx, dy=a.y-cy, dz=a.z-cz;
  const cY=Math.cos(rot.ay), sY=Math.sin(rot.ay);
  const rx=dx*cY-dz*sY, rz=dx*sY+dz*cY;
  dx=rx; dz=rz;
  const cX=Math.cos(rot.ax), sX=Math.sin(rot.ax);
  const ry=dy*cX-dz*sX, rz2=dy*sX+dz*cX;
  return { sx: cw/2+dx*scale, sy: ch/2-dy*scale, depth: rz2 };
}

// Split CA atoms into contiguous SS runs
interface SSRun { type: SS; atoms: Atom[] }

function splitRuns(ca: Atom[], ssMap: Map<number,SS>): SSRun[] {
  if (!ca.length) return [];
  const runs: SSRun[] = [];
  let curType = ssMap.get(ca[0].resSeq)||'loop';
  let curRun: Atom[] = [ca[0]];
  for (let i=1;i<ca.length;i++) {
    if (ca[i].resSeq - ca[i-1].resSeq > 1) {
      runs.push({type:curType,atoms:curRun}); curRun=[ca[i]]; curType=ssMap.get(ca[i].resSeq)||'loop'; continue;
    }
    const ss = ssMap.get(ca[i].resSeq)||'loop';
    if (ss===curType||(ss==='loop'&&curRun.length<3)||(curType==='loop'&&ss!=='loop'&&curRun.length<=1)) {
      if(ss!=='loop') curType=ss; curRun.push(ca[i]);
    } else { runs.push({type:curType,atoms:curRun}); curRun=[ca[i]]; curType=ss; }
  }
  runs.push({type:curType,atoms:curRun});
  return runs;
}

// ═══════════════════════════════════════════════════════════════════════════
//  COLOR HELPERS
// ═══════════════════════════════════════════════════════════════════════════

function lighten(hex: string, amt: number): string {
  if (hex.startsWith('rgb')) return hex;
  let c = hex.replace('#','');
  if (c.length===3) c=c[0]+c[0]+c[1]+c[1]+c[2]+c[2];
  if (c.length!==6) return hex;
  const r=Math.min(255,parseInt(c.substring(0,2),16)+amt);
  const g=Math.min(255,parseInt(c.substring(2,4),16)+amt);
  const b=Math.min(255,parseInt(c.substring(4,6),16)+amt);
  return `rgb(${r},${g},${b})`;
}

function darken(hex: string, amt: number): string {
  if (hex.startsWith('rgb')) return hex;
  let c = hex.replace('#','');
  if (c.length===3) c=c[0]+c[0]+c[1]+c[1]+c[2]+c[2];
  if (c.length!==6) return hex;
  const r=Math.max(0,parseInt(c.substring(0,2),16)-amt);
  const g=Math.max(0,parseInt(c.substring(2,4),16)-amt);
  const b=Math.max(0,parseInt(c.substring(4,6),16)-amt);
  return `rgb(${r},${g},${b})`;
}

function rgba(hex: string, a: number): string {
  let c = hex.replace('#','');
  if (c.length===3) c=c[0]+c[0]+c[1]+c[1]+c[2]+c[2];
  if (c.length!==6) return `rgba(128,128,128,${a})`;
  return `rgba(${parseInt(c.substring(0,2),16)},${parseInt(c.substring(2,4),16)},${parseInt(c.substring(4,6),16)},${a})`;
}

// ═══════════════════════════════════════════════════════════════════════════
//  SPLINE MATH
// ═══════════════════════════════════════════════════════════════════════════

interface Pt { x: number; y: number }

function spline(pts: Pt[], seg=6): Pt[] {
  if (pts.length<2) return [...pts];
  const out: Pt[] = [];
  for (let i=0;i<pts.length-1;i++) {
    const p0=pts[Math.max(0,i-1)], p1=pts[i], p2=pts[Math.min(pts.length-1,i+1)], p3=pts[Math.min(pts.length-1,i+2)];
    for (let t=0;t<seg;t++) {
      const tt=t/seg, t2=tt*tt, t3=t2*tt;
      out.push({
        x: .5*((2*p1.x)+(-p0.x+p2.x)*tt+(2*p0.x-5*p1.x+4*p2.x-p3.x)*t2+(-p0.x+3*p1.x-3*p2.x+p3.x)*t3),
        y: .5*((2*p1.y)+(-p0.y+p2.y)*tt+(2*p0.y-5*p1.y+4*p2.y-p3.y)*t2+(-p0.y+3*p1.y-3*p2.y+p3.y)*t3),
      });
    }
  }
  out.push(pts[pts.length-1]);
  return out;
}

function perpNormals(pts: Pt[]): {nx:number;ny:number}[] {
  return pts.map((_,i) => {
    let dx:number,dy:number;
    if (i===0) { dx=pts[1].x-pts[0].x; dy=pts[1].y-pts[0].y; }
    else if (i===pts.length-1) { dx=pts[i].x-pts[i-1].x; dy=pts[i].y-pts[i-1].y; }
    else { dx=pts[i+1].x-pts[i-1].x; dy=pts[i+1].y-pts[i-1].y; }
    const len=Math.sqrt(dx*dx+dy*dy)||1;
    return { nx:-dy/len, ny:dx/len };
  });
}

// ═══════════════════════════════════════════════════════════════════════════
//  CPK ATOM DATA
// ═══════════════════════════════════════════════════════════════════════════

const CPK: Record<string,string> = { C:'#909090',N:'#3050F8',O:'#FF0D0D',S:'#FFFF30',H:'#FFFFFF',P:'#FF8000',FE:'#E06633',ZN:'#7D80B0',MG:'#8AFF00',MN:'#9C7AC7' };
const VDW: Record<string,number> = { H:1.2,C:1.7,N:1.55,O:1.52,S:1.8,P:1.8 };

function atomColor(elem: string, cc: string) { return CPK[elem.toUpperCase()]||cc; }
function vdw(elem: string) { return VDW[elem.toUpperCase()]||1.7; }

// ═══════════════════════════════════════════════════════════════════════════
//  RENDERER
// ═══════════════════════════════════════════════════════════════════════════

function render(ctx: CanvasRenderingContext2D, model: Model, opts: {
  rot: Rot; hl: HighlightedRes[]; selChain?: string;
  mode: DisplayMode; labels: boolean; w: number; h: number; scale: number;
}) {
  const { rot, hl, selChain, mode, labels, w, h, scale } = opts;
  const { center, caByChain, chains, ssByChain, atoms, resByChain } = model;
  const [cx,cy,cz] = center;

  ctx.fillStyle = BG;
  ctx.fillRect(0,0,w,h);

  const hlSet = new Set<string>();
  for (const r of hl) hlSet.add(`${r.chain}:${r.res_num}`);

  const p = (a:Atom) => proj(a,cx,cy,cz,scale,rot,w,h);

  // ─── CARTOON ────────────────────────────────────────────────────────────
  if (mode === 'cartoon') {
    const segs: { type:SS; depth:number; pts:Pt[]; dim:boolean; cc:string; isHL:boolean }[] = [];

    for (const ch of chains) {
      const ca = caByChain.get(ch)||[];
      if (ca.length<2) continue;
      const dim = !!selChain && ch!==selChain;
      const ss = ssByChain.get(ch);
      const cc = chainColor(ch);
      const runs = ss ? splitRuns(ca,ss) : [{type:'loop' as SS,atoms:ca}];

      for (const run of runs) {
        if (run.atoms.length<2) continue;
        const pr = run.atoms.map(a=>p(a));
        const pts = pr.map(r=>({x:r.sx,y:r.sy}));
        const depth = pr.reduce((s,r)=>s+r.depth,0)/pr.length;
        const isHL = run.atoms.some(a=>hlSet.has(`${ch}:${a.resSeq}`));
        segs.push({type:run.type,depth,pts,dim,cc,isHL});
      }
    }

    segs.sort((a,b)=>a.depth-b.depth);

    for (const s of segs) {
      if (s.type==='helix') drawHelix(ctx,s);
      else if (s.type==='sheet') drawSheet(ctx,s);
      else drawLoop(ctx,s);
    }

    // Ligands
    for (const a of atoms) {
      if (!a.het) continue;
      if (selChain && a.chain!==selChain) continue;
      const r = p(a);
      ctx.save();
      ctx.globalAlpha=0.9;
      ctx.beginPath(); ctx.arc(r.sx,r.sy,4,0,Math.PI*2);
      ctx.fillStyle=LIGAND_CLR; ctx.fill();
      ctx.strokeStyle=darken(LIGAND_CLR,40); ctx.lineWidth=0.8; ctx.stroke();
      ctx.restore();
    }

    if (labels) for (const ch of chains) { if(selChain&&ch!==selChain) continue; drawLabels(ctx,caByChain.get(ch)||[],hlSet,ch,cx,cy,cz,scale,rot,w,h); }
    drawSSLegend(ctx,w,h);
    drawChainLegend(ctx,chains,h);
    return;
  }

  // ─── BALL & STICK ───────────────────────────────────────────────────────
  if (mode === 'ballstick') {
    const items: { t:'b'|'a'; depth:number; x1:number;y1:number;x2?:number;y2?:number;col:string;r?:number;hl:boolean }[] = [];

    const rMap = new Map<string,Map<number,Atom[]>>();
    for (const a of atoms) {
      if(a.het) continue;
      if(!rMap.has(a.chain)) rMap.set(a.chain,new Map());
      const m=rMap.get(a.chain)!;
      if(!m.has(a.resSeq)) m.set(a.resSeq,[]);
      m.get(a.resSeq)!.push(a);
    }

    for (const ch of chains) {
      const cm=rMap.get(ch); if(!cm) continue;
      const dim=!!selChain&&ch!==selChain;
      const cc=chainColor(ch);

      for (const [,ra] of cm) {
        for (let i=0;i<ra.length;i++) for (let j=i+1;j<ra.length;j++) {
          const d=dist(ra[i],ra[j]);
          if (ra[i].elem.toUpperCase()==='H'&&ra[j].elem.toUpperCase()==='H') continue;
          if(d>2||d<0.4) continue;
          const pa=p(ra[i]),pb=p(ra[j]);
          const isHL=hlSet.has(`${ch}:${ra[i].resSeq}`);
          items.push({t:'b',depth:(pa.depth+pb.depth)/2,x1:pa.sx,y1:pa.sy,x2:pb.sx,y2:pb.sy,col:dim?rgba(cc,0.2):isHL?HIGHLIGHT_CLR:'#808080',hl:isHL&&!dim});
        }
      }

      const ca=caByChain.get(ch)||[];
      for(let i=0;i<ca.length-1;i++){
        if(ca[i+1].resSeq-ca[i].resSeq!==1)continue;
        const pa=p(ca[i]),pb=p(ca[i+1]);
        const isHL=hlSet.has(`${ch}:${ca[i].resSeq}`)||hlSet.has(`${ch}:${ca[i+1].resSeq}`);
        items.push({t:'b',depth:(pa.depth+pb.depth)/2,x1:pa.sx,y1:pa.sy,x2:pb.sx,y2:pb.sy,col:dim?rgba(cc,0.2):isHL?HIGHLIGHT_CLR:'#808080',hl:isHL&&!dim});
      }

      for (const a of (resByChain.get(ch)||[])) {
        if(a.het) continue;
        const pr=p(a);
        const isHL=hlSet.has(`${ch}:${a.resSeq}`);
        items.push({t:'a',depth:pr.depth,x1:pr.sx,y1:pr.sy,col:dim?rgba(cc,0.2):isHL?HIGHLIGHT_CLR:atomColor(a.elem,cc),r:dim?1.5:vdw(a.elem)*1.6*Math.min(scale/8,1.2),hl:isHL&&!dim});
      }
    }

    for(const a of atoms){
      if(!a.het) continue;
      const pr=p(a);
      items.push({t:'a',depth:pr.depth,x1:pr.sx,y1:pr.sy,col:LIGAND_CLR,r:vdw(a.elem)*2*Math.min(scale/8,1.2),hl:true});
    }

    items.sort((a,b)=>a.depth-b.depth);

    for(const it of items){
      ctx.save();
      if(it.t==='b'&&it.x2!==undefined&&it.y2!==undefined){
        ctx.beginPath(); ctx.moveTo(it.x1,it.y1); ctx.lineTo(it.x2!,it.y2!);
        ctx.strokeStyle=it.col; ctx.lineWidth=2; ctx.lineCap='round';
        ctx.globalAlpha=it.hl?1:0.7; ctx.stroke();
      } else if(it.t==='a'){
        const r=Math.max(it.r||3,1.5);
        const g=ctx.createRadialGradient(it.x1-r*0.3,it.y1-r*0.3,r*0.05,it.x1,it.y1,r);
        g.addColorStop(0,lighten(it.col,60)); g.addColorStop(0.5,it.col); g.addColorStop(1,darken(it.col,40));
        ctx.beginPath(); ctx.arc(it.x1,it.y1,r,0,Math.PI*2);
        ctx.fillStyle=g;
        if(it.hl){ctx.shadowColor='rgba(0,0,0,0.25)';ctx.shadowBlur=3;}
        ctx.fill();
      }
      ctx.restore();
    }

    if(labels) for(const ch of chains){if(selChain&&ch!==selChain)continue;drawLabels(ctx,caByChain.get(ch)||[],hlSet,ch,cx,cy,cz,scale,rot,w,h);}
    drawChainLegend(ctx,chains,h);
    return;
  }

  // ─── BACKBONE ───────────────────────────────────────────────────────────
  const items: { t:'l'|'c'; depth:number; x1:number;y1:number;x2?:number;y2?:number;r?:number;col:string;hl:boolean }[] = [];

  for (const ch of chains) {
    const ca=caByChain.get(ch)||[];
    const cc=chainColor(ch);
    const dim=!!selChain&&ch!==selChain;

    for(let i=0;i<ca.length-1;i++){
      if(ca[i+1].resSeq-ca[i].resSeq>1)continue;
      const pa=p(ca[i]),pb=p(ca[i+1]);
      const isHL=hlSet.has(`${ch}:${ca[i].resSeq}`)||hlSet.has(`${ch}:${ca[i+1].resSeq}`);
      items.push({t:'l',depth:(pa.depth+pb.depth)/2,x1:pa.sx,y1:pa.sy,x2:pb.sx,y2:pb.sy,col:isHL?HIGHLIGHT_CLR:cc,hl:isHL&&!dim});
    }
    for(const a of ca){
      const pr=p(a);
      const isHL=hlSet.has(`${ch}:${a.resSeq}`);
      items.push({t:'c',depth:pr.depth,x1:pr.sx,y1:pr.sy,r:isHL?2.5:1.5,col:isHL?HIGHLIGHT_CLR:cc,hl:isHL&&!dim});
    }
    if(!dim) for(const a of atoms){ if(!a.het||a.chain!==ch)continue; const pr=p(a); items.push({t:'c',depth:pr.depth,x1:pr.sx,y1:pr.sy,r:4,col:LIGAND_CLR,hl:true}); }
  }

  items.sort((a,b)=>a.depth-b.depth);

  for(const it of items){
    ctx.save();
    if(it.t==='l'){
      ctx.beginPath(); ctx.moveTo(it.x1,it.y1); if(it.x2!==undefined&&it.y2!==undefined)ctx.lineTo(it.x2!,it.y2!);
      ctx.strokeStyle=it.col; ctx.lineWidth=3; ctx.lineCap='round';
      ctx.globalAlpha=it.hl?1:0.85;
      if(it.hl){ctx.shadowColor=it.col;ctx.shadowBlur=4;}
      ctx.stroke();
    } else if(it.t==='c'){
      const r=it.r||2;
      ctx.beginPath(); ctx.arc(it.x1,it.y1,r,0,Math.PI*2);
      ctx.fillStyle=it.col; ctx.globalAlpha=it.hl?1:0.8;
      if(it.hl){ctx.shadowColor=it.col;ctx.shadowBlur=4;}
      ctx.fill();
    }
    ctx.restore();
  }

  drawChainLegend(ctx,chains,h);
}

// ═══════════════════════════════════════════════════════════════════════════
//  CARTOON DRAWING HELPERS
// ═══════════════════════════════════════════════════════════════════════════

function drawHelix(ctx: CanvasRenderingContext2D, s: { pts:Pt[]; dim:boolean; cc:string; isHL:boolean }) {
  const {pts,dim,cc,isHL}=s;
  if(pts.length<2) return;
  const tw=dim?3.5:5.5;
  const base=dim?rgba(cc,0.25):HELIX_COLOR;
  const light=dim?rgba(cc,0.15):HELIX_LIGHT;
  const dark=dim?rgba(cc,0.35):HELIX_DARK;

  ctx.save(); ctx.globalAlpha=dim?0.35:1;

  const sm=spline(pts,5); if(sm.length<2){ctx.restore();return;}
  const ns=perpNormals(sm);
  const L:Pt[]=[], R:Pt[]=[];
  for(let i=0;i<sm.length;i++){
    L.push({x:sm[i].x+ns[i].nx*tw, y:sm[i].y+ns[i].ny*tw});
    R.push({x:sm[i].x-ns[i].nx*tw, y:sm[i].y-ns[i].ny*tw});
  }

  // Filled ribbon body
  ctx.beginPath(); ctx.moveTo(L[0].x,L[0].y);
  for(let i=1;i<L.length;i++) ctx.lineTo(L[i].x,L[i].y);
  for(let i=R.length-1;i>=0;i--) ctx.lineTo(R[i].x,R[i].y);
  ctx.closePath();

  // Gradient for 3D cylindrical look
  const mi=Math.floor(sm.length/2);
  const gx0=sm[mi].x+ns[mi].nx*tw, gy0=sm[mi].y+ns[mi].ny*tw;
  const gx1=sm[mi].x-ns[mi].nx*tw, gy1=sm[mi].y-ns[mi].ny*tw;

  if(isFinite(gx0)&&isFinite(gy0)&&isFinite(gx1)&&isFinite(gy1)&&(Math.abs(gx0-gx1)>0.5||Math.abs(gy0-gy1)>0.5)){
    const g=ctx.createLinearGradient(gx0,gy0,gx1,gy1);
    g.addColorStop(0,light); g.addColorStop(0.35,base); g.addColorStop(0.7,dark); g.addColorStop(1,darken(base,60));
    ctx.fillStyle=g;
  } else { ctx.fillStyle=base; }

  if(isHL&&!dim){ctx.shadowColor=HIGHLIGHT_CLR;ctx.shadowBlur=6;}
  ctx.fill();

  // Outline
  ctx.strokeStyle=darken(base,55); ctx.lineWidth=0.6; ctx.shadowBlur=0; ctx.globalAlpha=dim?0.12:0.2; ctx.stroke();

  // Specular highlight
  if(!dim){
    ctx.beginPath(); ctx.moveTo(L[0].x,L[0].y);
    for(let i=1;i<L.length;i++) ctx.lineTo(L[i].x,L[i].y);
    ctx.strokeStyle='rgba(255,255,255,0.35)'; ctx.lineWidth=1.2; ctx.globalAlpha=0.6; ctx.stroke();
  }

  // End caps
  for(const cap of [sm[0],sm[sm.length-1]]){
    ctx.beginPath(); ctx.arc(cap.x,cap.y,tw,0,Math.PI*2);
    ctx.fillStyle=base; ctx.globalAlpha=dim?0.35:1; ctx.fill();
  }

  ctx.restore();
}

function drawSheet(ctx: CanvasRenderingContext2D, s: { pts:Pt[]; dim:boolean; cc:string; isHL:boolean }) {
  const {pts,dim,cc,isHL}=s;
  if(pts.length<2) return;
  const aw=dim?5:8.5;
  const base=dim?rgba(cc,0.25):SHEET_COLOR;
  const light=dim?rgba(cc,0.15):SHEET_LIGHT;
  const dark=dim?rgba(cc,0.35):SHEET_DARK;

  ctx.save(); ctx.globalAlpha=dim?0.35:1;

  const sm=spline(pts,5); if(sm.length<3){ctx.restore();return;}
  const ns=perpNormals(sm);
  const bEnd=sm.length-Math.max(2,Math.ceil(14/2));
  if(bEnd<=1){ctx.restore();return;}

  const L:Pt[]=[], R:Pt[]=[];
  for(let i=0;i<bEnd;i++){
    L.push({x:sm[i].x+ns[i].nx*aw, y:sm[i].y+ns[i].ny*aw});
    R.push({x:sm[i].x-ns[i].nx*aw, y:sm[i].y-ns[i].ny*aw});
  }

  const tip=sm[sm.length-1];
  ctx.beginPath(); ctx.moveTo(L[0].x,L[0].y);
  for(let i=1;i<L.length;i++) ctx.lineTo(L[i].x,L[i].y);
  ctx.lineTo(tip.x,tip.y);
  for(let i=R.length-1;i>=0;i--) ctx.lineTo(R[i].x,R[i].y);
  ctx.closePath();

  const mi=Math.floor(bEnd/2);
  const gx0=sm[mi].x+ns[mi].nx*aw, gy0=sm[mi].y+ns[mi].ny*aw;
  const gx1=sm[mi].x-ns[mi].nx*aw, gy1=sm[mi].y-ns[mi].ny*aw;

  if(isFinite(gx0)&&isFinite(gy0)&&isFinite(gx1)&&isFinite(gy1)&&(Math.abs(gx0-gx1)>0.5||Math.abs(gy0-gy1)>0.5)){
    const g=ctx.createLinearGradient(gx0,gy0,gx1,gy1);
    g.addColorStop(0,light); g.addColorStop(0.45,base); g.addColorStop(1,dark);
    ctx.fillStyle=g;
  } else { ctx.fillStyle=base; }

  if(isHL&&!dim){ctx.shadowColor=HIGHLIGHT_CLR;ctx.shadowBlur=6;}
  ctx.fill();
  ctx.strokeStyle=darken(base,50); ctx.lineWidth=0.7; ctx.shadowBlur=0; ctx.globalAlpha=dim?0.12:0.25; ctx.stroke();

  if(!dim){
    ctx.beginPath(); ctx.moveTo(L[0].x,L[0].y);
    for(let i=1;i<L.length;i++) ctx.lineTo(L[i].x,L[i].y);
    ctx.lineTo(tip.x,tip.y);
    ctx.strokeStyle='rgba(255,255,255,0.3)'; ctx.lineWidth=1; ctx.globalAlpha=0.5; ctx.stroke();
  }

  ctx.restore();
}

function drawLoop(ctx: CanvasRenderingContext2D, s: { pts:Pt[]; dim:boolean; cc:string; isHL:boolean }) {
  const {pts,dim,cc,isHL}=s;
  if(pts.length<2) return;
  const tw=dim?1.5:2.5;
  const base=dim?rgba(cc,0.2):LOOP_COLOR;
  const light=dim?rgba(cc,0.1):LOOP_LIGHT;
  const dark=dim?rgba(cc,0.3):LOOP_DARK;

  ctx.save(); ctx.globalAlpha=dim?0.3:0.9;

  const sm=spline(pts,4); if(sm.length<2){ctx.restore();return;}
  const ns=perpNormals(sm);
  const L:Pt[]=[], R:Pt[]=[];
  for(let i=0;i<sm.length;i++){
    L.push({x:sm[i].x+ns[i].nx*tw, y:sm[i].y+ns[i].ny*tw});
    R.push({x:sm[i].x-ns[i].nx*tw, y:sm[i].y-ns[i].ny*tw});
  }

  ctx.beginPath(); ctx.moveTo(L[0].x,L[0].y);
  for(let i=1;i<L.length;i++) ctx.lineTo(L[i].x,L[i].y);
  for(let i=R.length-1;i>=0;i--) ctx.lineTo(R[i].x,R[i].y);
  ctx.closePath();

  const mi=Math.floor(sm.length/2);
  const gx0=sm[mi].x+ns[mi].nx*tw, gy0=sm[mi].y+ns[mi].ny*tw;
  const gx1=sm[mi].x-ns[mi].nx*tw, gy1=sm[mi].y-ns[mi].ny*tw;

  if(isFinite(gx0)&&isFinite(gy0)&&isFinite(gx1)&&isFinite(gy1)&&(Math.abs(gx0-gx1)>0.2||Math.abs(gy0-gy1)>0.2)){
    const g=ctx.createLinearGradient(gx0,gy0,gx1,gy1);
    g.addColorStop(0,light); g.addColorStop(0.4,base); g.addColorStop(1,dark);
    ctx.fillStyle=g;
  } else { ctx.fillStyle=base; }

  if(isHL&&!dim){ctx.shadowColor=HIGHLIGHT_CLR;ctx.shadowBlur=4;}
  ctx.fill();

  for(const cap of [sm[0],sm[sm.length-1]]){
    ctx.beginPath(); ctx.arc(cap.x,cap.y,tw,0,Math.PI*2);
    ctx.fillStyle=base; ctx.fill();
  }

  ctx.restore();
}

// ═══════════════════════════════════════════════════════════════════════════
//  CANVAS OVERLAY HELPERS
// ═══════════════════════════════════════════════════════════════════════════

function drawLabels(ctx: CanvasRenderingContext2D, ca:Atom[], hlSet:Set<string>, ch:string, cx:number,cy:number,cz:number,scale:number,rot:Rot,w:number,h:number) {
  ctx.save(); ctx.font='9px sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
  for(let i=0;i<ca.length;i++){
    if(i%5!==0) continue;
    const a=ca[i], pr=proj(a,cx,cy,cz,scale,rot,w,h);
    const txt=`${a.resName}${a.resSeq}`;
    const tw=ctx.measureText(txt).width+6, th=14;
    ctx.fillStyle='rgba(255,255,255,0.88)'; ctx.fillRect(pr.sx-tw/2,pr.sy-th/2-10,tw,th);
    ctx.strokeStyle='rgba(0,0,0,0.1)'; ctx.lineWidth=0.5; ctx.strokeRect(pr.sx-tw/2,pr.sy-th/2-10,tw,th);
    ctx.fillStyle='#374151'; ctx.fillText(txt,pr.sx,pr.sy-10);
  }
  ctx.restore();
}

function drawSSLegend(ctx: CanvasRenderingContext2D, w: number, h: number) {
  ctx.save(); ctx.font='10px sans-serif'; ctx.textBaseline='middle';
  const items=[{c:HELIX_COLOR,l:'α-Helix'},{c:SHEET_COLOR,l:'β-Sheet'},{c:LOOP_COLOR,l:'Loop'}];
  let lx=w-10; const ly=h-16;
  ctx.textAlign='right';
  for(let i=items.length-1;i>=0;i--){
    const it=items[i]; const tw=ctx.measureText(it.l).width;
    ctx.fillStyle='#6b7280'; ctx.fillText(it.l,lx,ly); lx-=tw+6;
    ctx.fillStyle=it.c; ctx.fillRect(lx-14,ly-4,12,8);
    ctx.strokeStyle=darken(it.c,30); ctx.lineWidth=0.5; ctx.strokeRect(lx-14,ly-4,12,8);
    lx-=22;
  }
  ctx.restore();
}

function drawChainLegend(ctx: CanvasRenderingContext2D, chains: string[], h: number) {
  ctx.save(); ctx.font='10px sans-serif'; ctx.textBaseline='middle';
  let lx=10; const ly=h-16;
  for(const ch of chains){
    ctx.fillStyle=chainColor(ch); ctx.beginPath(); ctx.arc(lx+4,ly,4,0,Math.PI*2); ctx.fill();
    ctx.fillStyle='#6b7280'; ctx.textAlign='left'; ctx.fillText(`Chain ${ch}`,lx+12,ly);
    lx+=65;
  }
  ctx.restore();
}

// ═══════════════════════════════════════════════════════════════════════════
//  REACT COMPONENT
// ═══════════════════════════════════════════════════════════════════════════

export default function PdbViewer3D({
  pdbText, highlightedResidues=[], selectedChain, height='650px', className, onResidueClick,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const modelRef = useRef<Model | null>(null);
  const rafRef = useRef<number>(0);
  const [modelInfo, setModelInfo] = useState<{ atoms: number; chains: number } | null>(null);

  const [mode, setMode] = useState<DisplayMode>('cartoon');
  const [spin, setSpin] = useState(false);
  const [labels, setLabels] = useState(false);
  const [rot, setRot] = useState<Rot>({ ax: -0.4, ay: 0.6 });
  const [zoom, setZoom] = useState(1);
  const [dragging, setDragging] = useState(false);
  const [lastMouse, setLastMouse] = useState({ x: 0, y: 0 });

  // Parse PDB text
  useEffect(() => {
    if (!pdbText) { modelRef.current = null; React.startTransition(() => setModelInfo(null)); return; }
    const m = parsePdb(pdbText);
    modelRef.current = m;
    React.startTransition(() => setModelInfo({ atoms: m.atoms.length, chains: m.chains.length }));
  }, [pdbText]);

  // Auto-fit zoom
  useEffect(() => {
    const m = modelRef.current;
    if (!m || !containerRef.current) return;
    const rect = containerRef.current.getBoundingClientRect();
    const fitScale = Math.min(rect.width, rect.height) / (m.size * 1.1) * 0.8;
    React.startTransition(() => setZoom(fitScale > 0 ? fitScale : 1));
  }, [pdbText]);

  // Render loop
  useEffect(() => {
    const canvas = canvasRef.current;
    const m = modelRef.current;
    if (!canvas || !m) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    let running = true;

    const frame = () => {
      if (!running) return;

      // Spin
      if (spin) {
        setRot(prev => ({ ax: prev.ax, ay: prev.ay + 0.008 }));
      }

      const rect = containerRef.current?.getBoundingClientRect();
      if (rect) {
        canvas.width = rect.width * window.devicePixelRatio;
        canvas.height = rect.height * window.devicePixelRatio;
        canvas.style.width = rect.width + 'px';
        canvas.style.height = rect.height + 'px';
        ctx.setTransform(window.devicePixelRatio, 0, 0, window.devicePixelRatio, 0, 0);
      }

      render(ctx, m, {
        rot, hl: highlightedResidues, selChain: selectedChain,
        mode, labels, w: rect?.width||600, h: rect?.height||650, scale: zoom,
      });

      rafRef.current = requestAnimationFrame(frame);
    };

    frame();
    return () => { running = false; cancelAnimationFrame(rafRef.current); };
  }, [rot, zoom, mode, labels, spin, highlightedResidues, selectedChain, pdbText]);

  // Mouse handlers
  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    setDragging(true);
    setLastMouse({ x: e.clientX, y: e.clientY });
  }, []);

  const handleMouseMove = useCallback((e: React.MouseEvent) => {
    if (!dragging) return;
    const dx = e.clientX - lastMouse.x;
    const dy = e.clientY - lastMouse.y;
    setRot(prev => ({ ax: prev.ax + dy * 0.005, ay: prev.ay + dx * 0.005 }));
    setLastMouse({ x: e.clientX, y: e.clientY });
  }, [dragging, lastMouse]);

  const handleMouseUp = useCallback(() => {
    setDragging(false);
  }, []);

  const handleWheel = useCallback((e: React.WheelEvent) => {
    e.preventDefault();
    setZoom(prev => Math.max(0.1, Math.min(50, prev * (1 - e.deltaY * 0.001))));
  }, []);

  // Residue click (backbone mode)
  const handleCanvasClick = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!onResidueClick || !modelRef.current) return;
    const canvas = canvasRef.current;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const m = modelRef.current;
    const [cx,cy,cz] = m.center;

    let closest: { chain:string; resSeq:number; resName:string; d:number } | null = null;
    for (const ch of m.chains) {
      for (const a of (m.caByChain.get(ch)||[])) {
        const pr = proj(a,cx,cy,cz,zoom,rot,rect.width,rect.height);
        const d = Math.sqrt((pr.sx-mx)**2+(pr.sy-my)**2);
        if (d<15 && (!closest || d<closest.d)) closest={chain:ch,resSeq:a.resSeq,resName:a.resName,d};
      }
    }
    if (closest) onResidueClick({chain:closest.chain, res_num:closest.resSeq, res_name:closest.resName});
  }, [onResidueClick, zoom, rot]);

  return (
    <div ref={containerRef} className={`relative ${className||''}`} style={{height}}>
      {/* Canvas */}
      <canvas
        ref={canvasRef}
        className="absolute inset-0 cursor-grab active:cursor-grabbing"
        style={{width:'100%',height:'100%'}}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onWheel={handleWheel}
        onClick={handleCanvasClick}
      />

      {/* Toolbar */}
      <div className="absolute top-2 left-2 flex flex-wrap gap-1 z-10">
        {MODES.map(m => (
          <button
            key={m.id}
            onClick={() => setMode(m.id)}
            className={`px-2 py-1 rounded text-[10px] font-medium transition-all border ${
              mode===m.id
                ? 'bg-emerald-600 text-white border-emerald-600 shadow-sm'
                : 'bg-white/90 text-gray-600 border-gray-200 hover:border-emerald-300 hover:bg-emerald-50'
            }`}
          >
            {m.icon} {m.label}
          </button>
        ))}
        <button
          onClick={() => setSpin(s=>!s)}
          className={`px-2 py-1 rounded text-[10px] font-medium transition-all border ${
            spin
              ? 'bg-violet-600 text-white border-violet-600 shadow-sm'
              : 'bg-white/90 text-gray-600 border-gray-200 hover:border-violet-300 hover:bg-violet-50'
          }`}
        >
          {spin ? '⏸ Spin On' : '▶ Spin'}
        </button>
        <button
          onClick={() => setLabels(l=>!l)}
          className={`px-2 py-1 rounded text-[10px] font-medium transition-all border ${
            labels
              ? 'bg-amber-500 text-white border-amber-500 shadow-sm'
              : 'bg-white/90 text-gray-600 border-gray-200 hover:border-amber-300 hover:bg-amber-50'
          }`}
        >
          {labels ? '🏷 Labels On' : '🏷 Labels'}
        </button>
        <button
          onClick={() => { setRot({ax:-0.4,ay:0.6}); const m=modelRef.current; if(m&&containerRef.current){const r=containerRef.current.getBoundingClientRect(); setZoom(Math.min(r.width,r.height)/(m.size*1.1)*0.8);} }}
          className="px-2 py-1 rounded text-[10px] font-medium bg-white/90 text-gray-600 border border-gray-200 hover:border-gray-300 hover:bg-gray-50 transition-all"
        >
          ↺ Reset
        </button>
      </div>

      {/* Info badge */}
      {modelInfo && (
        <div className="absolute bottom-2 right-2 px-2 py-1 bg-white/80 backdrop-blur-sm rounded text-[9px] text-gray-400 border border-gray-100">
          {modelInfo.atoms} atoms · {modelInfo.chains} chain{modelInfo.chains>1?'s':''} · Drag to rotate · Scroll to zoom
        </div>
      )}
    </div>
  );
}
