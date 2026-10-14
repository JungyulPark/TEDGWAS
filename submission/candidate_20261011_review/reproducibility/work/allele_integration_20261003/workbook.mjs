import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {FileBlob,SpreadsheetFile} from '@oai/artifact-tool';
process.on('uncaughtException',e=>{console.error(e.message);process.exit(1);});
const W=path.dirname(fileURLToPath(import.meta.url));
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(path.join(W,'edit_subset.xlsx')));
const preview=async(sheetName,range,name)=>{
 const b=await wb.render({sheetName,range,scale:1,format:'png'});
 await fs.writeFile(path.join(W,name+'.png'),new Uint8Array(await b.arrayBuffer()));
};
if(process.argv.includes('--before')){
 console.log((await wb.inspect({kind:'table',range:'Readme!A1:B6',include:'values',tableMaxRows:6,tableMaxCols:2,maxChars:1800})).ndjson);
 await preview('Readme','A1:B8','workbook_before');
 process.exit(0);
}
const blocks=JSON.parse(await fs.readFile(path.join(W,'blocks.json'),'utf8'));
const rm=wb.worksheets.getItem('Readme'),dict=wb.worksheets.getItem('Column_dictionary');
rm.getRange('B1').values=[['Five-outcome TSHR/IGF1R reappraisal with allele audit — 2026-10-03']];
rm.getRange('A25:B25').values=[['TSHR_allele_audit','Post hoc audit: blood alleles, signed reference LD, five disease datasets, MR reproduction and all inferred haplotypes. No new tissue MR.']];
rm.getRange('A25:B25').format.font={name:'Calibri',size:11,color:'#000000'};
rm.getRange('A25:B25').format.wrapText=true;rm.getRange('A25:B25').format.rowHeight=36;
const sh=wb.worksheets.add('TSHR_allele_audit');sh.showGridLines=false;
sh.getRange('A1:AA78').format.font={name:'Arial',size:10,color:'#000000'};
sh.getRange('A1:AA78').format.columnWidth=19;
sh.getRange('A1:AA78').format.rowHeight=20;
sh.getRange('A1').values=[['TSHR allele direction audit']];sh.getRange('A1').format.font={size:14,bold:true,color:'#000000'};
sh.getRange('A2').values=[['Sections A–E preserve analytical exports. Alleles, scales and source P values remain explicit.']];
sh.getRange('A3').values=[['Two linked variants are not independent replications. Haplotype frequencies are inferred from unphased reference genotypes.']];
const definitions={
 'Pvalue':'Original blood eQTL P value; Zscore refers to AssessedAllele.',
 'AssessedAllele':'Allele to which the signed blood-expression Z statistic refers.',
 'OtherAllele':'Alternative allele at the same variant.',
 'Zscore':'Signed association Z statistic for blood expression.',
 'NrSamples':'Variant-specific eQTLGen expression sample size.',
 'SNPPos':'GRCh37 SNP coordinate, not interchangeable with disease GRCh38 coordinates.',
 'GenePos':'Source eQTLGen GRCh37 gene-position anchor.',
 'signed_G_dosage_r':'Pearson correlation of G-allele copy counts at both named variants.',
 'r_squared':'Squared genotype-count correlation, not the haplotype maximum-likelihood r-squared.',
 'phase_basis':'Whether coupling is inferred rather than directly observed phased haplotypes.',
 'source_beta':'Disease coefficient in the source file before sign alignment and, for UKB, scale conversion.',
 'source_se':'Standard error in the original source scale.',
 'source_p':'Original disease GWAS P value, not overwritten by a normal-approximation MR P.',
 'source_effect_scale':'UKB uses linear-model coefficients; BBJ and FinnGen use log odds.',
 'aligned_G_beta_logOR':'Disease log-odds coefficient aligned to G after the specified scale conversion.',
 'aligned_G_se':'Disease standard error on the aligned log-odds scale.',
 'aligned_G_OR':'Per-G-allele disease odds ratio, not the expression-scaled MR odds ratio.',
 'aligned_G_CI_low':'Lower 95% confidence limit for per-G-allele disease OR.',
 'aligned_G_CI_high':'Upper 95% confidence limit for per-G-allele disease OR.',
 'aligned_G_eaf':'Frequency of G, after aligning the original effect-allele frequency.',
 'sign_multiplier':'1 for retained source direction; -1 for the opposite source effect allele.',
 'recomputed_normal_p':'Two-sided normal-approximation P computed from aligned disease beta/SE.',
 'MR_OR':'Reproduced OR per approximately one SD of reconstructed blood expression.',
 'MR_p':'Reproduced normal-approximation Wald P, distinct from source GWAS P.',
 'CI_z_multiplier':'Historical CI convention: 1.96 or qnorm(0.975), preserved per original pipeline.',
 'source_gwas_p':'Original GWAS P reported separately from the reconstructed MR P.',
 'exposure_F':'Instrument strength calculated as the squared original expression Z statistic.',
 'max_numeric_difference':'Largest absolute difference from canonical beta, SE, OR or CI.',
 'estimated_frequency':'PLINK two-locus maximum-likelihood frequency, rounded in its original log output.',
 'haplotype':'rs179252 allele followed by rs179247 allele; coupling inferred from unphased genotypes.',
 'plink_raw_A1_correlation':'PLINK correlation using its A1 allele at each variant; compare with the explicit G/G correlation.',
 'n_missing':'Reference participants missing at least one of the two genotypes.',
};
const dd=[];
for(const b of blocks){
 const n=b.headers.length,r=b.row;
 sh.getRangeByIndexes(r-1,0,1,1).values=[[b.title]];
 sh.getRangeByIndexes(r-1,0,1,1).format.font={bold:true,color:'#000000'};
 const head=sh.getRangeByIndexes(r,0,1,n);head.values=[b.headers];
 head.format={fill:'#E8EEF4',font:{bold:true,color:'#000000'},wrapText:true,rowHeight:60,verticalAlignment:'center',horizontalAlignment:'center'};
 const values=b.values.map(row=>row.map((v,c)=>b.headers[c]==='method'?'PLINK ML, unphased':b.headers[c]==='phase_basis'?'Inferred from unphased genotypes':v));
 const body=sh.getRangeByIndexes(r+1,0,b.values.length,n);body.values=values;
 body.format.verticalAlignment='center';body.setNumberFormat('0.000000');
 for(let c=0;c<n;c++){
  const h=b.headers[c],col=sh.getRangeByIndexes(r+1,c,b.values.length,1);
  if(/(^Pvalue$|BonferroniP|^FDR$|^source_p$|^MR_p$|recomputed_normal_p|source_gwas_p|difference)/.test(h))col.setNumberFormat('0.000E+00');
  else if(/^(n_|Nr|.*_pos37|.*pos38|SNPPos|GenePos|SNPChr|GeneChr)|^sign_multiplier$/.test(h))col.setNumberFormat('0');
  else if(b.values.every(v=>v[c]===null||typeof v[c]!=='number'))col.setNumberFormat('General');
  if(h==='outcome'||h==='method'||h==='phase_basis')col.format.wrapText=true;
  dd.push(['TSHR_allele_audit',`${b.title[0]}: ${h}`,definitions[h]||`Original exported field in section ${b.title}. See the section-specific source and Supplementary Methods.`]);
 }
 body.format.rowHeight=b.title.startsWith('B ')?48:32;
 sh.getRangeByIndexes(r+b.values.length+2,0,1,1).values=[['Source: '+b.file]];
}
dict.getRangeByIndexes(278,0,dd.length,3).values=dd;
dict.getRangeByIndexes(278,0,dd.length,3).format.font={name:'Calibri',size:11,color:'#000000'};
dict.getRangeByIndexes(278,0,dd.length,3).format.wrapText=true;
dict.getRangeByIndexes(278,0,dd.length,3).format.rowHeight=42;
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'TSHR_allele_audit!A5:J8',include:'values,formulas',tableMaxRows:4,tableMaxCols:10,maxChars:2400})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!',options:{useRegex:true,maxResults:5},maxChars:800})).ndjson);
await preview('Readme','A23:B25','readme_updated');
await preview('Column_dictionary','A279:C284','dictionary_updated');
for(const b of blocks){await preview('TSHR_allele_audit',`A${b.row}:J${b.row+Math.min(b.values.length+2,7)}`,'audit_'+b.title[0]);
 if(b.headers.length>10)await preview('TSHR_allele_audit',`K${b.row+1}:${String.fromCharCode(64+b.headers.length)}${b.row+Math.min(b.values.length+1,6)}`,'audit_'+b.title[0]+'_right');
}
await (await SpreadsheetFile.exportXlsx(wb)).save(path.join(W,'edited_subset.xlsx'));
console.log(JSON.stringify({auditBlocks:blocks.length,newDictionaryRows:dd.length,exported:true}));
