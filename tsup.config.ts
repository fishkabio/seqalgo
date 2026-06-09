import { defineConfig } from 'tsup';

export default defineConfig({
  entry: {
    index: 'index.ts',
    'needle/index': 'src/needle/index.ts',
  },
  format: ['esm', 'cjs'],
  dts: true,
  minify: true,
  sourcemap: true,
  clean: true,
  splitting: false,
  treeshake: true,
  outDir: 'dist',
});
