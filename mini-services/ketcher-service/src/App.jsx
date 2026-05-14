import React, { useEffect, useRef, useCallback } from 'react';
// Import from binaryWasm entry point - uses file-based Worker URLs
// (default entry uses WorkerFactory blob URLs which break WASM loading)
import { StandaloneStructServiceProvider } from 'ketcher-standalone/dist/binaryWasm';
import 'ketcher-react/dist/index.css';

let ketcherInstance = null;

function App() {
  const containerRef = useRef(null);
  const ketcherRef = useRef(null);

  useEffect(() => {
    // Listen for messages from parent iframe
    const handleMessage = async (event) => {
      const { type, id, data } = event.data || {};

      switch (type) {
        case 'getSmiles': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'getSmilesResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            const smiles = await ketcherRef.current.getSmiles();
            window.parent.postMessage({ type: 'getSmilesResult', id, data: smiles }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'getSmilesResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'getMolfile': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'getMolfileResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            const molfile = await ketcherRef.current.getMolfile();
            window.parent.postMessage({ type: 'getMolfileResult', id, data: molfile }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'getMolfileResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'getKet': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'getKetResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            const ket = await ketcherRef.current.getKet();
            window.parent.postMessage({ type: 'getKetResult', id, data: ket }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'getKetResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'setMolecule': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'setMoleculeResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            await ketcherRef.current.setMolecule(data);
            window.parent.postMessage({ type: 'setMoleculeResult', id, success: true }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'setMoleculeResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'clear': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'clearResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            await ketcherRef.current.setMolecule('');
            window.parent.postMessage({ type: 'clearResult', id, success: true }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'clearResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'getRxn': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'getRxnResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            const rxn = await ketcherRef.current.getRxn();
            window.parent.postMessage({ type: 'getRxnResult', id, data: rxn }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'getRxnResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'automap':
        case 'automapCustom': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'automapResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            // Get current structure as KET, automap via Indigo, set back in editor
            const ket = await ketcherRef.current.getKet();
            const structService = ketcherRef.current.structService;
            const fmt = (data && data.output_format) || 'chemical/x-indigo-ket';
            const result = await structService.automap({
              struct: ket,
              output_format: fmt,
              mode: 'alter',
            });
            // Update the editor with the automapped structure
            await ketcherRef.current.setMolecule(result.struct);
            window.parent.postMessage({ type: 'automapResult', id, data: result.struct }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'automapResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'automapDirect': {
          // Bypass Ketcher structService, use Indigo directly
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'automapResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            const ket = await ketcherRef.current.getKet();
            const structService = ketcherRef.current.structService;
            // Try with empty output_format (Indigo may default to something useful)
            const result = await structService.automap({
              struct: ket,
              mode: 'alter',
            });
            await ketcherRef.current.setMolecule(result.struct);
            // Get SMILES after automap
            const smiles = await ketcherRef.current.getSmiles();
            window.parent.postMessage({ type: 'automapResult', id, data: smiles }, '*');
          } catch (e) {
            window.parent.postMessage({ type: 'automapResult', id, error: e.message }, '*');
          }
          break;
        }
        case 'generateImage': {
          if (!ketcherRef.current) {
            window.parent.postMessage({ type: 'generateImageResult', id, error: 'Ketcher not initialized' }, '*');
            return;
          }
          try {
            const blob = await ketcherRef.current.generateImage('mol', {
              width: data?.width || 500,
              height: data?.height || 400,
            });
            const reader = new FileReader();
            reader.onload = () => {
              window.parent.postMessage({ type: 'generateImageResult', id, data: reader.result }, '*');
            };
            reader.readAsDataURL(blob);
          } catch (e) {
            window.parent.postMessage({ type: 'generateImageResult', id, error: e.message }, '*');
          }
          break;
        }
      }
    };

    window.addEventListener('message', handleMessage);
    return () => window.removeEventListener('message', handleMessage);
  }, []);

  const handleInit = useCallback((ketcher) => {
    ketcherRef.current = ketcher;
    ketcherInstance = ketcher;
    // Notify parent that Ketcher is ready
    window.parent.postMessage({ type: 'ketcherReady' }, '*');
  }, []);

  const handleError = useCallback((message) => {
    console.error('Ketcher error:', message);
  }, []);

  const structServiceProvider = useRef(new StandaloneStructServiceProvider());

  // Dynamically import and render the Editor component
  const [EditorComponent, setEditorComponent] = React.useState(null);

  useEffect(() => {
    import('ketcher-react').then((mod) => {
      setEditorComponent(() => mod.Editor);
    }).catch((err) => {
      console.error('Failed to load ketcher-react:', err);
    });
  }, []);

  if (!EditorComponent) {
    return (
      <div style={{
        width: '100%',
        height: '100%',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        fontFamily: 'system-ui, sans-serif',
        color: '#666',
      }}>
        <div style={{ textAlign: 'center' }}>
          <div style={{
            width: 40,
            height: 40,
            border: '3px solid #e0e0e0',
            borderTopColor: '#00bfa5',
            borderRadius: '50%',
            animation: 'spin 1s linear infinite',
            margin: '0 auto 12px'
          }} />
          <style>{`@keyframes spin { to { transform: rotate(360deg); } }`}</style>
          Loading Chemical Editor...
        </div>
      </div>
    );
  }

  return (
    <div ref={containerRef} style={{ width: '100%', height: '100%' }}>
      <EditorComponent
        staticResourcesUrl={new URL('.', window.location.href).pathname.replace(/\/$/, '') || '/'}
        structServiceProvider={structServiceProvider.current}
        onInit={handleInit}
        errorHandler={handleError}
        ketcherId="editor-1"
      />
    </div>
  );
}

export default App;
