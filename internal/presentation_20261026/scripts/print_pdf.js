// node print_pdf.js page.html out.pdf  (A4, footer with page numbers)
const { chromium } = require('playwright');
const path = require('path');
(async () => {
  const [src, out] = process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto('file://' + path.resolve(src), { waitUntil: 'networkidle' });
  await page.pdf({
    path: out, format: 'A4', printBackground: true, displayHeaderFooter: true,
    headerTemplate: '<div></div>',
    footerTemplate: '<div style="font-family:Noto Sans KR,sans-serif;font-size:7pt;color:#5B6472;width:100%;padding:0 15mm;' +
      'display:flex;justify-content:space-between;"><span>발표자료 수정 보고서 · 2026-10-26</span>' +
      '<span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>',
    margin: { top: '16mm', bottom: '16mm', left: '15mm', right: '15mm' },
  });
  await browser.close();
})();
