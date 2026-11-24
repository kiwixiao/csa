#!/usr/bin/env node
// Word version of the CFD report: one US Letter page per report page image (the same PNGs the
// PDF pages are made from, 8.5 x 11 in at 300 dpi), each scaled to fit the page margins.
//   node report_to_docx.js --title "ENT005 CFD report" out.docx page1.png page2.png ...
// needs the npm package 'docx' (set NODE_PATH if it is not installed next to this script)
const fs = require("fs");
const { Document, Packer, Paragraph, ImageRun, AlignmentType } = require("docx");

const args = process.argv.slice(2);
let title = "CFD report";
if (args[0] === "--title") { title = args[1]; args.splice(0, 2); }
const [out, ...pngs] = args;
if (!out || pngs.length === 0) { console.error("usage: report_to_docx.js [--title T] out.docx page.png ..."); process.exit(1); }

const IN = 1440;                                   // DXA per inch
const PAGE_W = 8.5, PAGE_H = 11, MARGIN = 0.3;     // inches
const PX = 96;                                     // docx-js image size unit: px at 96 dpi

function pngSize(buf) { return { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) }; }

const pages = pngs.map((p, i) => {
  const data = fs.readFileSync(p);
  const { w, h } = pngSize(data);
  const boxW = PAGE_W - 2 * MARGIN, boxH = PAGE_H - 2 * MARGIN - 0.15;   // a little slack for the line height
  const s = Math.min(boxW / w, boxH / h);
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    pageBreakBefore: i > 0,
    spacing: { before: 0, after: 0 },
    children: [new ImageRun({ type: "png", data, transformation: { width: Math.round(w * s * PX), height: Math.round(h * s * PX) },
                              altText: { title: `${title}, page ${i + 1}`, description: `${title}, page ${i + 1}`, name: `page${i + 1}` } })],
  });
});

const doc = new Document({
  title,
  sections: [{
    properties: { page: { size: { width: PAGE_W * IN, height: PAGE_H * IN },
                          margin: { top: MARGIN * IN, bottom: MARGIN * IN, left: MARGIN * IN, right: MARGIN * IN } } },
    children: pages,
  }],
});
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(out, buf); console.log(`wrote ${out} (${pngs.length} pages)`); });
