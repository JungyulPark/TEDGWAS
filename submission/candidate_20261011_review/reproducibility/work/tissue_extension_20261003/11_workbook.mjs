import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {FileBlob,SpreadsheetFile} from '@oai/artifact-tool';
process.on('uncaughtException',e=>{console.error(e.message);process.exit(1);});
const W=path.dirname(fileURLToPath(import.meta.url));
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(path.join(W,'edit_subset.xlsx')));
const blocks=JSON.parse(await fs.readFile(path.join(W,'blocks.json'),'utf8'));
const rm=wb.worksheets.getItem('Readme'),dict=wb.worksheets.getItem('Column_dictionary');
rm.getRange('B1').values=[['Five-outcome reappraisal with allele audit and fixed GTEx tissue extension — 2026-10-03']];
rm.getRange('A26:B26').values=[['Tissue_extension','All 50 planned tissue MR comparisons, ten-combination availability, four selected instruments and complete candidate/alignment records. Table S3 now presents tissue results; original detection benchmarks remain in Detection_original.']];
rm.getRange('A26:B26').format={font:{name:'Calibri',size:11,color:'#000000'},wrapText:true,rowHeight:64};
const sh=wb.worksheets.add('Tissue_extension');sh.showGridLines=false;
sh.getRange('A1:AJ153').format={font:{name:'Arial',size:10,color:'#000000'},columnWidth:19,rowHeight:24};
sh.getRange('A1').values=[['TSHR and IGF1R tissue extension']];sh.getRange('A1').format.font={bold:true,size:14,color:'#000000'};
sh.getRange('A2').values=[['Frozen panel: two genes in five GTEx v8 tissues. Ten primary and forty secondary comparisons; unavailable cells retained.']];
sh.getRange('A3').values=[['Statistical-program exports. One instrument per estimable combination; no tissue colocalization. These cells are not editable model inputs.']];
const defs={P:'Two-sided normal Wald MR P value; not the source eQTL or GWAS P.',P_adjusted:'min(1, family size × P): 10 for primary FinnGen GD, 40 for secondary disease comparisons.',status:'estimated or not_estimable; a missing estimate is not evidence of no effect.',reason:'Specific cause of an unestimable planned comparison.',OR:'Disease odds ratio per one unit of normalized inverse-normal-transformed expression; do not compare effect magnitudes directly to reconstructed eQTLGen units.',CI_lower:'Lower 95% Wald confidence limit for the MR OR.',CI_upper:'Upper 95% Wald confidence limit for the MR OR.',direction_vs_blood:'Comparison with known blood direction, TSHR lower odds and IGF1R higher odds; direction alone does not establish association.',minimum_eQTL_P:'Smallest nominal eQTL P in the all-tested-gene eGenes table; blank when not tested.',tested:'Whether the gene is present in the tissue eGenes file containing all tested genes.',n_exposure:'GTEx v8 tissue QTL sample count from official metadata, checked against allele count and MAF for selected variants.',exposure_eaf:'Actual ALT frequency from selected eGene MAF and ref_factor; never MAF without allele orientation.',pval_nominal:'Original GTEx nominal eQTL P; used for the fixed P < 5e-8 selection threshold.',slope:'GTEx original ALT-allele effect on normalized expression.',slope_se:'Original standard error of GTEx slope.',F:'Squared slope / slope_se ratio.',in_LD_reference:'Whether exact rsID and compatible GRCh37 position/alleles were available in the EUR clumping panel.',genome_build:'Coordinate build for position fields; pos38 is GRCh38 and explicitly named pos37 is GRCh37.',Q:'Unavailable: every estimable tissue MR uses one SNP; no heterogeneity test.',Q_df:'Unavailable for one-SNP Wald MR.',Q_pvalue:'Unavailable for one-SNP Wald MR.',coloc_status:'Not run because complete regional summaries were not obtained; significant-only eQTLs were never used.',reference_effect_allele_frequency:'EUR reference frequency, distinct from the GTEx tissue effect-allele frequency.',source_p_outcome:'Original GWAS P preserved separately from the Wald MR P.',scale_factor:'One-time UKB conversion 1/[mu(1-mu)], mu=3731/484598; one for BBJ/FinnGen.',allele_sign:'Sign multiplier aligning source disease effect allele to the GTEx ALT allele.',availability_status:'Candidate availability before LD clumping; n_selected_instruments gives final retained exposure variants.'};
let dd=[];
for(const b of blocks){
 const n=b.headers.length,r=b.row;sh.getRangeByIndexes(r-1,0,1,1).values=[[b.title]];sh.getRangeByIndexes(r-1,0,1,1).format.font={bold:true,color:'#000000'};
 const head=sh.getRangeByIndexes(r,0,1,n);head.values=[b.headers];head.format={fill:'#E8EEF4',font:{bold:true,color:'#000000'},wrapText:true,rowHeight:62,verticalAlignment:'center',horizontalAlignment:'center'};
 const body=sh.getRangeByIndexes(r+1,0,b.values.length,n);body.values=b.values;body.format={wrapText:true,rowHeight:50,verticalAlignment:'center'};
 for(let j=0;j<n;j++){
  const k=b.headers[j],range=sh.getRangeByIndexes(r+1,j,b.values.length,1);
  const numeric=b.values.every(v=>v[j]===null||typeof v[j]==='number');
  if(numeric)range.setNumberFormat(/(^P$|^P_adjusted$|_P$|pval|^p_|qvalue|source_p)/i.test(k)?'0.000E+00':/^(n_|pos|chromosome|ma_count|ma_samples|tss_distance|planned_family_tests|regional_variants|allele_sign)/.test(k)?'0':'0.000000');
  else range.setNumberFormat('General');
  dd.push(['Tissue_extension',b.title[0]+': '+k,defs[k]||`Exported ${k} in section ${b.title[0]}; source ${b.file}. See fixed-panel Supplementary Methods.`]);
 }
 sh.getRangeByIndexes(r+b.values.length+2,0,1,1).values=[['Source: '+b.file]];
}
dict.getRangeByIndexes(353,0,dd.length,3).values=dd;dict.getRangeByIndexes(353,0,dd.length,3).format={font:{name:'Calibri',size:11,color:'#000000'},wrapText:true,rowHeight:44};
wb.recalculate();console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!',options:{useRegex:true,maxResults:5},maxChars:600})).ndjson);
for(const [range,name] of [['A5:J12','availability_left'],['K6:P12','availability_right'],['A19:L27','results_top'],['G28:N36','results_lower'],['O20:AG25','results_detail'],['A73:J78','instruments_left'],['K74:AC78','instruments_right'],['A81:J87','harmonized_left'],['K82:AB87','harmonized_right'],['A105:J111','candidates_left'],['K106:AA111','candidates_right']]){
 const b=await wb.render({sheetName:'Tissue_extension',range,scale:1,format:'png'});await fs.writeFile(path.join(W,name+'.png'),new Uint8Array(await b.arrayBuffer()));
}
for(const [sheetName,range,name] of [['Readme','A25:B26','readme_updated'],['Column_dictionary','A354:C360','dictionary_updated']]){
 const b=await wb.render({sheetName,range,scale:1,format:'png'});await fs.writeFile(path.join(W,name+'.png'),new Uint8Array(await b.arrayBuffer()));
}
await (await SpreadsheetFile.exportXlsx(wb)).save(path.join(W,'edited_subset.xlsx'));console.log(JSON.stringify({dictionaryRows:dd.length,newSheet:'Tissue_extension'}));
