import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// Served as a subpage at https://ezrakruger.cc/originmarker/ -> relative base, never an absolute origin.
export default defineConfig({
  plugins: [react()],
  base: './',
  server: {
    proxy: {
      '/api': { target: 'http://127.0.0.1:8000', changeOrigin: true },
    },
  },
})
