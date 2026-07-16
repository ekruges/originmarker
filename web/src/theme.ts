import { createTheme, type MantineColorsTuple } from '@mantine/core'

// Ensembl/UCSC-ish: restrained blue, grey borders, dense type, near-square corners.
const genomeBlue: MantineColorsTuple = [
  '#eef4fa', '#dbe7f3', '#b7cee6', '#8fb2d7', '#6d9acb',
  '#578bc4', '#4a83c1', '#337ab7', '#2e6da4', '#255a88',
]

const genomeGrey: MantineColorsTuple = [
  '#f7f8f9', '#eceef0', '#dcdfe3', '#c6cbd1', '#adb4bc',
  '#8d959e', '#6b727b', '#4d545c', '#333940', '#1d2126',
]

export const theme = createTheme({
  primaryColor: 'genomeBlue',
  primaryShade: 8,
  colors: { genomeBlue, genomeGrey },
  white: '#ffffff',
  fontFamily:
    '-apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif',
  fontFamilyMonospace: 'ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace',
  defaultRadius: 2,
  radius: { xs: '1px', sm: '2px', md: '2px', lg: '3px', xl: '4px' },
  fontSizes: { xs: '11px', sm: '12px', md: '13px', lg: '14px', xl: '16px' },
  lineHeights: { xs: '1.3', sm: '1.35', md: '1.4', lg: '1.45', xl: '1.5' },
  headings: {
    fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif',
    sizes: {
      h1: { fontSize: '19px', fontWeight: '600', lineHeight: '1.3' },
      h2: { fontSize: '16px', fontWeight: '600', lineHeight: '1.3' },
      h3: { fontSize: '14px', fontWeight: '600', lineHeight: '1.3' },
      h4: { fontSize: '13px', fontWeight: '600', lineHeight: '1.3' },
    },
  },
  shadows: { xs: 'none', sm: 'none', md: 'none', lg: 'none', xl: 'none' },
  components: {
    Paper: { defaultProps: { withBorder: true, shadow: undefined, radius: 2 } },
    Button: { defaultProps: { radius: 2, size: 'xs' } },
    TextInput: { defaultProps: { radius: 2, size: 'xs' } },
    NumberInput: { defaultProps: { radius: 2, size: 'xs' } },
    Select: { defaultProps: { radius: 2, size: 'xs' } },
    Autocomplete: { defaultProps: { radius: 2, size: 'xs' } },
    // Mantine uppercases badge labels by default: caps read as an alarm, and these are
    // information.
    Badge: {
      defaultProps: { radius: 2 },
      styles: { label: { textTransform: 'none' as const, fontWeight: 600 } },
    },
    Alert: { defaultProps: { radius: 2 } },
    Table: { defaultProps: { horizontalSpacing: 6, verticalSpacing: 3, fontSize: '12px' } },
  },
})
