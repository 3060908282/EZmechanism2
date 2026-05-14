/**
 * Type declarations for 3Dmol.js loaded via script tag.
 *
 * The 3Dmol library is loaded as an external script (/3Dmol.min.js)
 * and exposes a global `$3Dmol` on the window object.
 * We declare its types here for TypeScript support.
 */

interface ViewerSpec {
  backgroundColor?: string;
  backgroundAlpha?: number;
  antialias?: boolean;
  disableFog?: boolean;
  nomouse?: boolean;
  hoverDuration?: number;
  cartoonQuality?: number;
  minimumZoomToDistance?: number;
  lowerZoomLimit?: number;
  upperZoomLimit?: number;
  callback?: (viewer: $3Dmol) => void;
}

interface StyleSpec {
  cartoon?: {
    color?: string;
    opacity?: number;
    style?: string;
    tubes?: boolean;
    arrows?: boolean;
    colorscheme?: string | Record<string, string>;
  };
  stick?: {
    radius?: number;
    color?: string;
    opacity?: number;
    colorscheme?: string | Record<string, string>;
    hidden?: boolean;
  };
  sphere?: {
    scale?: number;
    color?: string;
    opacity?: number;
    colorscheme?: string | Record<string, string>;
    hidden?: boolean;
  };
  line?: {
    color?: string;
    linewidth?: number;
    hidden?: boolean;
  };
  cross?: {
    linewidth?: number;
    color?: string;
    hidden?: boolean;
    radius?: number;
  };
  ballAndStick?: {
    color?: string;
    opacity?: number;
    colorscheme?: string | Record<string, string>;
  };
}

interface ModelSpec {
  hetflag?: boolean;
  chain?: string;
  resi?: number | [number, number] | string;
  resn?: string;
  atom?: string;
  elem?: string;
  model?: number;
  invert?: boolean;
  not?: ModelSpec;
  or?: ModelSpec[];
  byres?: boolean;
  expand?: number;
  within?: { distance: number; sel: ModelSpec };
  withinGroup?: { distance: number; sel: ModelSpec };
}

interface AtomSpec {
  resn?: string;
  resi?: number;
  chain?: string;
  atom?: string;
  elem?: string;
  x?: number;
  y?: number;
  z?: number;
  hetflag?: boolean;
  serial?: number;
  index?: number;
  ss?: string;
  color?: string;
  b?: number;
  bonds?: number[];
  bondOrder?: number[];
  clickable?: boolean;
  callback?: (atom: AtomSpec, viewer: $3Dmol, event: Event, container: HTMLElement) => void;
  hoverable?: boolean;
  hover_callback?: (atom: AtomSpec, viewer: $3Dmol, event: Event, container: HTMLElement) => void;
  unhover_callback?: (atom: AtomSpec, viewer: $3Dmol) => void;
  /** Dynamically attached label reference (used by hover system) */
  label?: any;
  hoverLabel?: any;
  [key: string]: any;
}

interface LabelSpec {
  position: { x: number; y: number; z: number };
  fontSize?: number;
  fontColor?: string;
  fontWeight?: string;
  backgroundColor?: string;
  backgroundOpacity?: number;
  borderColor?: string;
  borderOpacity?: number;
  borderThickness?: number;
  font?: string;
  inFront?: boolean;
  showBackground?: boolean;
  alignment?: string;
  borderRadius?: number;
}

interface SphereSpec {
  center: { x: number; y: number; z: number };
  radius: number;
  color?: string;
  opacity?: number;
}

class $3Dmol {
  static createViewer(
    element: HTMLElement | string,
    config?: ViewerSpec
  ): $3Dmol;

  addModel(data?: string, format?: string, options?: Record<string, unknown>): $3Dmol;
  addModels(data: string, format: string, options?: Record<string, unknown>): $3Dmol[];
  removeModel(model?: $3Dmol | number): $3Dmol;
  removeAllModels(): $3Dmol;
  removeAllSurfaces(): $3Dmol;
  removeAllShapes(): $3Dmol;
  removeAllLabels(): $3Dmol;

  setStyle(sel: ModelSpec | StyleSpec, style?: StyleSpec): $3Dmol;
  addStyle(sel: ModelSpec | StyleSpec, style?: StyleSpec): $3Dmol;

  zoomTo(sel?: ModelSpec, animationDuration?: number): $3Dmol;
  zoom(factor?: number, animationDuration?: number): $3Dmol;
  render(callback?: () => void): $3Dmol;
  resize(): $3Dmol;
  rotate(angle: number, axis?: string, animationDuration?: number): $3Dmol;
  center(sel?: ModelSpec, animationDuration?: number): $3Dmol;

  addLabel(text: string, options: LabelSpec): any;
  removeLabel(label: any): $3Dmol;
  addSphere(spec: SphereSpec): any;
  addArrow(spec: any): any;
  addCylinder(spec: any): any;
  addLine(spec: any): any;

  setClickable(
    sel: ModelSpec,
    clickable: boolean,
    callback: (atom: AtomSpec, viewer: $3Dmol, event: Event, container: HTMLElement, ...rest: any[]) => void
  ): $3Dmol;
  setHoverable(
    sel: ModelSpec,
    hoverable: boolean,
    hover_callback: (atom: AtomSpec, viewer: $3Dmol, event: Event, container: HTMLElement, ...rest: any[]) => void,
    unhover_callback: (atom: AtomSpec, viewer: $3Dmol) => void
  ): $3Dmol;

  enableContextMenu(sel: ModelSpec, enabled: boolean): $3Dmol;

  setColorByElement(sel: ModelSpec, colors: Record<string, string>): $3Dmol;
  setColorByProperty(sel: ModelSpec, prop: string, scheme: any): $3Dmol;

  selectedAtoms(sel: ModelSpec): AtomSpec[];
  pdbData(sel: ModelSpec): string;
  getUniqueValues(attribute: string, sel?: ModelSpec): string[];

  setBackgroundColor(color: string, alpha?: number): $3Dmol;
  setViewStyle(parameters: any): $3Dmol;
  setProjection(proj: string): $3Dmol;

  setHoverDuration(duration?: number): $3Dmol;
  spin(axis?: string | boolean, speed?: number): void;

  pngURI(): string;
  getCanvas(): HTMLCanvasElement;

  getModel(id?: number): $3Dmol;

  static SurfaceType: {
    VDW: number;
    MS: number;
    SAS: number;
    SES: number;
  };

  static ColorSpec: any;
  static Gradient: {
    RWB: any;
    Sinebow: any;
    Roygb: any;
  };
}

// Global declaration for script-tag loaded 3Dmol
interface Window {
  $3Dmol: typeof $3Dmol;
}
