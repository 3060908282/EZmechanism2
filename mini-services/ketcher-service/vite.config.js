import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  server: { port: 3004, host: '0.0.0.0' },
  base: '/ketcher/',
  build: { outDir: 'dist' },
  resolve: {
    alias: {
      // lodash uses "global" which doesn't exist in browsers
      global: 'global-jsx/polyfill',
    },
  },
  define: {
    // Provide global polyfill for browser
    global: 'globalThis',
  },
});
