const express = require('express');
const path = require('path');

const app = express();
const PORT = 3004;
const DIST_DIR = path.join(__dirname, 'dist');

// Serve static files from dist at the root
app.use(express.static(DIST_DIR));

// Serve ketcher-standalone WASM and worker files from node_modules
// The StandaloneStructServiceProvider needs to find these at the staticResourcesUrl
const ketcherStandaloneDir = path.join(__dirname, 'node_modules', 'ketcher-standalone', 'dist');
app.use('/binaryWasm', express.static(path.join(ketcherStandaloneDir, 'binaryWasm')));
app.use('/binaryWasmNoRender', express.static(path.join(ketcherStandaloneDir, 'binaryWasmNoRender')));

// Serve WASM files with correct MIME type
app.use('/assets', (req, res, next) => {
  if (req.path.endsWith('.wasm')) {
    res.type('application/wasm');
  }
  express.static(path.join(DIST_DIR, 'assets'))(req, res, next);
});

// Serve root - return index.html for any unmatched path (SPA behavior)
app.get('/', (req, res) => {
  res.sendFile(path.join(DIST_DIR, 'index.html'));
});

app.listen(PORT, '0.0.0.0', () => {
  console.log(`Ketcher editor service running on port ${PORT}`);
});
