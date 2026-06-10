import { defineConfig } from 'tsup';

export default defineConfig({
  entry: {
    index: 'index.ts',
    'needle/index': 'src/needle/index.ts',
    'sequence/index': 'src/sequence/index.ts',
    'alignment/index': 'src/alignment/index.ts',
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
