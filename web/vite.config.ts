import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// Served at the root on originmarker.app AND as a subpage at ezrakruger.cc/originmarker,
// from one build -> relative base, never an absolute origin.
export default defineConfig({
  plugins: [react()],
  base: './',
  server: {
    proxy: {
      '/api': { target: 'http://127.0.0.1:8000', changeOrigin: true },
    },
  },
})
