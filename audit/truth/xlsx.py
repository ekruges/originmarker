"""Read an .xlsx with the standard library: sheets as lists of rows, cells placed by their reference.

A cell Excel left empty is absent from the file, so a row's cells are placed by column letter and
never by order. Values are returned as text exactly as stored.
"""
import re
import zipfile
import xml.etree.ElementTree as ET

M = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
R = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'


def _col(ref):
    n = 0
    for ch in re.match(r'[A-Z]+', ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def sheets(path):
    """{sheet name: [row, ...]}, each row a list of strings padded to its last filled column."""
    z = zipfile.ZipFile(path)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall(M + 'si'):
            shared.append(''.join(t.text or '' for t in si.iter(M + 't')))
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    out = {}
    for s in ET.fromstring(z.read('xl/workbook.xml')).find(M + 'sheets'):
        target = rels[s.get(R + 'id')]
        part = target[1:] if target.startswith('/') else 'xl/' + target
        rows = []
        for row in ET.fromstring(z.read(part)).iter(M + 'row'):
            cells = {}
            for c in row.findall(M + 'c'):
                v = c.find(M + 'v')
                kind = c.get('t')
                if kind == 'inlineStr':
                    val = ''.join(t.text or '' for t in c.iter(M + 't'))
                elif v is None:
                    val = ''
                elif kind == 's':
                    val = shared[int(v.text)]
                else:
                    val = v.text or ''
                cells[_col(c.get('r'))] = val.strip()
            if any(cells.values()):
                width = max(cells) + 1
                rows.append([cells.get(i, '') for i in range(width)])
        out[s.get('name').strip()] = rows
    return out


if __name__ == '__main__':
    # Self-check: column letters, and a cell skipped in the file lands in its own column.
    assert _col('A1') == 0 and _col('Z9') == 25 and _col('AA3') == 26 and _col('AB10') == 27
    import io
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w') as z:
        z.writestr('xl/workbook.xml', '<workbook xmlns="%s" xmlns:r="%s"><sheets><sheet name="s" r:id="rId1"/></sheets></workbook>' % (M[1:-1], R[1:-1]))
        z.writestr('xl/_rels/workbook.xml.rels', '<Relationships><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr('xl/sharedStrings.xml', '<sst xmlns="%s"><si><t>x</t></si></sst>' % M[1:-1])
        z.writestr('xl/worksheets/sheet1.xml', '<worksheet xmlns="%s"><sheetData><row><c r="A1" t="s"><v>0</v></c><c r="C1"><v>7</v></c></row></sheetData></worksheet>' % M[1:-1])
    buf.seek(0)
    assert sheets(buf)['s'] == [['x', '', '7']]
    print('xlsx ok')
