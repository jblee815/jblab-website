(function(){
/* ---------- hero: DNA mismatch repair cycle on a double helix ---------- */
const cv=document.getElementById('field'),ctx=cv.getContext('2d');
let paused=false;
const stEls=[...document.querySelectorAll('#stages li')];
const COL={muts:'18,140,86',mutl:'204,118,14',exo:'200,46,62',pol:'32,104,196'};
const DUR=[4.5,1.6,6,3,3,1.2,1.8]; // search, recognition, signaling, excision, resynthesis, ligation, hold
const STAGE=[0,0,1,2,3,3,3];
const AMP=15,K=2*Math.PI/96,STEP=12;
let W,H,dpr,y0,xm,xn,xe,phase=0,t=0,ph=0,last=performance.now();

const ease=p=>p<.5?2*p*p:1-Math.pow(-2*p+2,2)/2;
const cl=p=>Math.max(0,Math.min(1,p));
function newCycle(){xm=Math.round((W*.68)/STEP)*STEP+6; // mismatch stays at the same position every cycle
const gap=Math.min(300,W*.4);xn=xm-gap;xe=xm+30;phase=0;t=0;}
function layout(){dpr=Math.min(devicePixelRatio||1,2);W=cv.clientWidth;H=cv.clientHeight;cv.width=W*dpr;cv.height=H*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);y0=H-82;
  if(xm===undefined){newCycle();return;}
  // on resize keep the running cycle; only move the fixed sites with the new width
  xm=Math.round((W*.68)/STEP)*STEP+6;xn=xm-Math.min(300,W*.4);xe=xm+30;}
const yT=x=>y0+AMP*Math.sin(K*x+ph), yB=x=>y0-AMP*Math.sin(K*x+ph);
function strand(f,x1,x2,col,w){if(x2<=x1)return;ctx.beginPath();for(let x=x1;x<=x2;x+=3)ctx[x===x1?'moveTo':'lineTo'](x,f(x));ctx.lineTo(x2,f(x2));ctx.strokeStyle=col;ctx.lineWidth=w;ctx.stroke();}
// ring(x, color, alpha, rx, ry, label, lineWidth, dy = vertical offset, below = label under the DNA)
function ring(x,c,a,rx,ry,label,lw,dy,below){if(a<=0)return;dy=dy||0;const y=y0+dy;ctx.save();ctx.globalAlpha=a;const g=ctx.createRadialGradient(x,y,0,x,y,ry*1.6);g.addColorStop(0,`rgba(${c},.28)`);g.addColorStop(1,`rgba(${c},0)`);ctx.fillStyle=g;ctx.fillRect(x-ry*1.6,y-ry*1.6,ry*3.2,ry*3.2);ctx.beginPath();ctx.ellipse(x,y,rx,ry,0,0,Math.PI*2);ctx.strokeStyle=`rgba(${c},.95)`;ctx.lineWidth=lw||3;ctx.shadowColor=`rgba(${c},.8)`;ctx.shadowBlur=4;ctx.stroke();ctx.shadowBlur=0;ctx.font='500 11px "IBM Plex Mono", monospace';ctx.fillStyle=`rgba(${c},1)`;ctx.textAlign='center';ctx.fillText(label,x,below?y+ry+16:y-ry-8);ctx.restore();}
// ATP binding flash: an expanding ring with a small ATP tag
function atp(x,c,k){if(k<=0||k>=1)return;ctx.save();ctx.globalAlpha=1-k;ctx.beginPath();ctx.arc(x,y0,8+30*k,0,Math.PI*2);ctx.strokeStyle=`rgba(${c},1)`;ctx.lineWidth=2;ctx.stroke();ctx.font='600 10px "IBM Plex Mono", monospace';ctx.fillStyle='rgba(20,32,29,1)';ctx.textAlign='center';ctx.fillText('+ATP',x,y0+48);ctx.restore();}
function draw(){
  ctx.fillStyle='#eaf1ee';ctx.fillRect(0,0,W,H);
  const p=cl(t/DUR[phase]);
  // gap in the daughter (bottom) strand
  let gL=0,gR=0;
  if(phase===3){gL=xn;gR=xn+(xe-xn)*ease(p);}
  if(phase===4){gL=xn+(xe-xn)*ease(p);gR=xe;}
  // original strand-discrimination nick at xn until excision starts; after resynthesis the only nick is at xe, where the new DNA meets the downstream strand, until ligase seals it
  const nickAt=(phase===2||phase===3)?xn:(phase===5&&p<.6)?xe:null;
  const repaired=phase>4||(phase===4&&xn+(xe-xn)*ease(p)>xm+4);
  // base-pair rungs
  for(let x=6;x<W;x+=STEP){
    if(gR>gL&&x>gL&&x<gR)continue;
    const mm=Math.abs(x-xm)<1;
    if(mm&&!repaired){const ym=(yT(x)+yB(x))/2;ctx.lineWidth=3;ctx.strokeStyle='rgba(200,46,62,1)';ctx.beginPath();ctx.moveTo(x-2,yT(x));ctx.lineTo(x-2,ym-3);ctx.stroke();ctx.strokeStyle='rgba(204,118,14,1)';ctx.beginPath();ctx.moveTo(x+2,yB(x));ctx.lineTo(x+2,ym+3);ctx.stroke();continue;}
    ctx.lineWidth=mm?3:2;ctx.strokeStyle=mm?'rgba(18,140,86,.95)':'rgba(60,90,80,.3)';ctx.beginPath();ctx.moveTo(x,yT(x));ctx.lineTo(x,yB(x));ctx.stroke();
  }
  // strands: template (top) intact, daughter (bottom) with nick and gap
  strand(yT,0,W,'rgba(38,56,51,.85)',2.4);
  const segs=[];let a=0;const cuts=[];
  if(nickAt!==null)cuts.push([nickAt-3,nickAt+3]);if(gR>gL)cuts.push([gL,gR]);cuts.sort((u,v)=>u[0]-v[0]);
  cuts.forEach(c=>{segs.push([a,Math.min(c[0],W)]);a=Math.max(a,c[1]);});segs.push([a,W]);
  segs.forEach(s=>strand(yB,s[0],s[1],'rgba(38,56,51,.85)',2.4));
  // new DNA grows continuously from the upstream 3′ end at xn (no break there)
  if(phase===4){strand(yB,xn,gL,'rgba(32,104,196,.9)',2.4);}
  if(phase>4){strand(yB,xn,nickAt!==null?xe-3:xe,`rgba(32,104,196,${phase===5?.9:Math.max(0,.9-.9*p)})`,2.4);}
  // mismatch glow
  if(!repaired){const g=ctx.createRadialGradient(xm,y0,0,xm,y0,34);g.addColorStop(0,'rgba(200,46,62,.35)');g.addColorStop(1,'rgba(200,46,62,0)');ctx.fillStyle=g;ctx.fillRect(xm-34,y0-34,68,68);}
  else if(phase>=5){const al=phase===5?.5:.5*(1-p);const g=ctx.createRadialGradient(xm,y0,0,xm,y0,34);g.addColorStop(0,`rgba(18,140,86,${al})`);g.addColorStop(1,'rgba(18,140,86,0)');ctx.fillStyle=g;ctx.fillRect(xm-34,y0-34,68,68);}
  // proteins: faint trail of recent random-walk positions
  const tt=performance.now()/1000;
  if(phase===0){const xs=Math.min(40,xn);const x=xs+(xm-xs)*ease(p)+26*Math.sin(tt*6.3)*(1-p)+14*Math.sin(tt*11)*(1-p);ring(x,COL.muts,1,17,32,'MutS');}
  if(phase===1){ring(xm,COL.muts,1,17+2*Math.sin(tt*9),32+2*Math.sin(tt*9),'MutS');}
  if(phase===2){
    // (a) MutS binds ATP at the mismatch and becomes a sliding clamp
    const a=cl(p/.18),sw=ease(a);
    // (b) the MutS·ATP clamp leaves the mismatch toward the nick
    const q=cl((p-.18)/.82),e=ease(q);
    const xs=xm+(xn+26-xm)*e+(q>0?6*Math.sin(tt*8)*(1-q):0);
    ring(xs,COL.muts,1,17-4*sw,32-6*sw,a<.35?'MutS':'MutS·ATP',3+.5*sw);
    atp(xm,COL.muts,a<1?a:0);
    // (c) once the MutS·ATP clamp is well clear of the mismatch, MutL comes straight down
    //     from solution and docks directly onto MutS·ATP (it never touches the mismatch site)
    const bind=cl((q-.30)/.20),xl=xs-14; // docks on the side facing away from the mismatch
    const la=cl((q-.54)/.14);
    if(q>=.30&&la<=0)ring(xl,COL.mutl,bind,13,26,'MutL',3.5,-70*(1-ease(bind)),bind>=.95);
    // (d) only after docking does MutL bind ATP, forming MutL·ATP on the moving complex
    if(la>0){atp(xl,COL.mutl,la<1?la:0);ring(xl,COL.mutl,1,13-1.5*ease(la),26-2*ease(la),la<.4?'MutL':'MutL·ATP',3.5,0,true);}
    ring(xn-16,COL.exo,cl((q-.72)*3.5),12,22,'Exonuclease',3);
  }
  if(phase===3){const fx=1-cl(p*1.6);ring(xn+26,COL.muts,fx,13,26,'MutS·ATP',3.5);ring(xn+12,COL.mutl,fx,11.5,24,'MutL·ATP',3.5,0,true);ring(xn+(xe-xn)*ease(p),COL.exo,1,12,22,'Exonuclease',3);}
  if(phase===4){ring(xn+(xe-xn)*ease(p),COL.pol,1,14,24,'Pol δ',3);}
  if(phase===5){const r=6+30*p;ctx.beginPath();ctx.arc(xe,y0,r,0,Math.PI*2);ctx.strokeStyle=`rgba(20,32,29,${1-p})`;ctx.lineWidth=2;ctx.stroke();ctx.font='500 11px "IBM Plex Mono", monospace';ctx.fillStyle=`rgba(20,32,29,${1-p*.6})`;ctx.textAlign='center';ctx.fillText('Ligase',xe,y0-34);}
  // 5′ / 3′ strand ends
  ctx.font='11px "IBM Plex Mono", monospace';ctx.fillStyle='rgba(83,101,95,.9)';ctx.textAlign='left';ctx.fillText('3′',4,y0-AMP-10);ctx.fillText('5′',4,y0+AMP+20);ctx.textAlign='right';ctx.fillText('5′',W-4,y0-AMP-10);ctx.fillText('3′',W-4,y0+AMP+20);
  stEls.forEach((li,i)=>li.classList.toggle('on',i===STAGE[phase]));
}
function tick(now){const dt=Math.min(.05,(now-last)/1000);last=now;if(paused||!W){requestAnimationFrame(tick);return;}t+=dt;
  if(t>=DUR[phase]){t=0;phase++;if(phase>=DUR.length)newCycle();}
  draw();requestAnimationFrame(tick);}
layout();draw();requestAnimationFrame(tick);
// redraw whenever the hero changes size (covers a zero-size canvas at first load)
let lastW=0,lastH=0;const relayout=()=>{if(cv.clientWidth&&(cv.clientWidth!==lastW||cv.clientHeight!==lastH)){lastW=cv.clientWidth;lastH=cv.clientHeight;layout();draw();}};
if('ResizeObserver' in window)new ResizeObserver(relayout).observe(cv);
addEventListener('resize',relayout);addEventListener('load',relayout);
const btn=document.getElementById('anim-btn');
btn.addEventListener('click',()=>{paused=!paused;btn.textContent=paused?'Play animation':'Pause animation';btn.setAttribute('aria-pressed',paused);});

/* ---------- research ---------- */
const T={
dna:[
 {k:'MutS · MutL · MMR',h:'DNA mismatch repair',p:'MutS homologs search for and recognize a mismatch, the signal is passed downstream with MLH/PMS and ATP, the error-containing strand is excised from a 5′ or 3′ end, and the gap is resynthesized and ligated. Defects in human MSH and MLH genes cause Lynch syndrome. Real-time single-molecule imaging lets us resolve the mechanics of each step.',r:[],img:['img/mmr-model.webp','Model of E. coli mismatch repair: MutS searches the DNA as an ADP-bound clamp by helical diffusion, recognizes a mismatch, exchanges ADP for ATP and becomes a sliding clamp that diffuses fast without ATP hydrolysis. MutL loads onto MutS, and MutH joins at the hemimethylated site.',2436,1050],side:true}
],
cell:[
 {k:'TNT · DFB',h:'Tunneling nanotubes',p:'Tunneling nanotubes are sub-micrometer membrane routes that directly connect distant cells and carry organelles, vesicles, RNAs and receptors. We proposed how close-ended TNTs form: HILO imaging reveals double filopodial bridges, STORM and SRRF show they twist into a double helix, and torsional and shearing energy drive the bridge-to-nanotube transition.',r:[['Chang et al., Science Advances 8, eabj3995 (2022)','https://www.science.org/doi/10.1126/sciadv.abj3995']],img:['img/tnt-model.webp','Model of close-ended TNT formation: filopodia from two cells meet through weak N-cadherin binding (DFB0), grow and twist into a double filopodial bridge attached to the paired cell (DFB1), then N-cadherin bonds break and the released filopodium retracts, leaving one filopodium held by N-cadherin clusters as a close-ended TNT.',1554,1170]}
],
tech:[
 {k:'Skybridge',h:'DNA Skybridge',p:'DNA molecules are suspended between the tops of ~4 µm-high quartz posts, so the surface below stays out of focus. Light totally reflected inside the quartz forms a thin sheet along the suspended DNA, exciting only the proteins bound to it, such as Cy5-labeled PCNA. Hundreds of DNA molecules can be imaged at once with little background, for high-throughput single-molecule studies of DNA–protein interactions.',r:[['Kim et al., Nucleic Acids Research 47, e107 (2019)','https://doi.org/10.1093/nar/gkz625']],
  media:[{type:'img',src:'img/skybridge.webp',w:828,h:1218,alt:'Top: electron micrograph of quartz posts 4.0 µm high and 13.1 µm apart with DNA drawn suspended between them. Bottom: schematic of DNA tethered across two posts, with a light sheet at the top exciting Cy5-PCNA on the DNA while the quartz bottom stays dark.'},
         {type:'video',src:'img/skybridge-movie.mp4',poster:'img/skybridge-poster.jpg',w:512,h:512,cap:'Movie of Cy5-PCNA on lambda phage DNA suspended on DNA skybridges',legend:[['green','Green: Lambda phage DNA'],['red','Red: Cy5-PCNA']]}]},
 {k:'Deep learning',h:'smDeepFLUOR',p:'smDeepFLUOR classifies the binding state of a molecule from how a single fluorophore looks and flickers. Each detected spot becomes a 7 \u00d7 7 pixel, 10-frame input unit that carries both its spatial intensity distribution and its temporal fluctuation. A 3D convolutional network extracts features from these units and a dense layer assigns the binding state. Trained on mNG-PABP spots, it separated specific from nonspecific binding with 93% and 97% accuracy on held-out test data, even though the two intensity distributions overlap heavily.',r:[['Lee et al., Nature Communications 17, 7810 (2026)','https://www.nature.com/articles/s41467-026-74716-3']],wide:true,
  media:[{type:'img',src:'img/smdeepfluor-pipeline.webp',w:2000,h:1542,alt:'smDeepFLUOR workflow: single-fluorophore imaging produces 7 by 7 pixel, 10-frame input units carrying spatial distribution and temporal fluctuation; three 3D convolution layers with max-pooling extract features, and a dense layer of 384, 16 and 2 neurons classifies them.'},
         {type:'img',src:'img/smdeepfluor-result.webp',w:1189,h:889,alt:'Classification of mNG-PABP binding: the full dataset is split 80% training and 20% test; intensity distributions of specific and nonspecific binding overlap, yet the confusion matrix shows 0.97 correct for nonspecific and 0.93 correct for specific binding.'}]}
]};
const topicsEl=document.getElementById('topics');
// lab members (current and alumni), including the spellings used in author lists
const LAB=['Jong-Bong Lee','Byungju Kim','Minsun Choi','Songmi Oh','Philsang Yu','Phil Sang Yu','ManSeong Yoon','Keunsang Yang','GyeongHo Kang',
 'Jinseob Lee','Inho Yang','Jincheol Seol','Gayun Bu','Jaeho Oh','Daehyung Kim','Ryanggeun Lee','Minhyeok Chang','Minhyoek Chang','Yongmoon Jeon',
 'Jong-Hyun Park','Jonghyun Park','Seung Hwan Lee','Sung-Sil Lee','Seongsil Lee','Won-Ki Cho','Soo-Ran Kim','Sooran Kim','Dong-Wook Kim','Daekil In','Jungsic Oh','Cherlhyun Jeong','C. Jeong'];
// a name matches with or without its hyphen/space and with an en dash, e.g. Jong–Bong Lee, Won Ki Cho
const labRe=new RegExp('(?<![\\w-])('+[...new Set(LAB.map(n=>n.split(/[- ]/)).map(p=>p.slice(0,-1).join('[-– ]?')+' '+p[p.length-1]))].sort((a,b)=>b.length-a.length).join('|')+')([*#]{0,2})(?![\\w-])','gi');
function auth(a){return esc(a).replace(labRe,'<strong>$1$2</strong>')}
function esc(s){return s.replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));}
function showTopics(g){
  document.querySelectorAll('[data-g]').forEach(b=>b.setAttribute('aria-selected',b.dataset.g===g));
  // portrait figures (taller than wide-ish) sit beside the text; landscape figures span the card below it
  topicsEl.innerHTML=T[g].map(t=>{
    let fig=t.img?`<figure class="topic-fig"><img src="${t.img[0]}" alt="${esc(t.img[1])}" width="${t.img[2]}" height="${t.img[3]}"></figure>`:'';
    if(t.media)fig=`<div class="topic-media">${t.media.map(m=>m.type==='video'
      ?`<figure class="topic-fig vid"><video src="${m.src}" poster="${m.poster}" width="${m.w}" height="${m.h}" autoplay muted loop playsinline controls preload="metadata" aria-label="${esc(m.cap)}"></video><figcaption class="legend">${(m.legend||[]).map(([c,l])=>`<span><i class="dot-${c}"></i>${esc(l)}</span>`).join('')}</figcaption></figure>`
      :`<figure class="topic-fig"><img src="${m.src}" alt="${esc(m.alt)}" width="${m.w}" height="${m.h}"></figure>`).join('')}</div>`;
    const refs=t.r.length?`<div class="refs">${t.r.map(r=>`<a href="${r[1]}">${esc(r[0])} ↗</a>`).join('')}</div>`:'';
    const text=`<div class="topic-text"><span class="k">${esc(t.k)}</span><h3>${esc(t.h)}</h3><p>${esc(t.p)}</p>${refs}</div>`;
    const side=!t.stack&&(t.side||t.media||(t.img&&t.img[2]/t.img[3]<1.6));
    return side?`<article class="topic side${t.wide?' wide':''}">${text}${fig}</article>`:`<article class="topic">${text}${fig}</article>`;
  }).join('');
}
document.querySelectorAll('[data-g]').forEach(b=>b.addEventListener('click',()=>showTopics(b.dataset.g)));
showTopics('dna');

/* ---------- people ---------- */
const init=n=>n.split(/[\s-]+/).filter(Boolean).map(w=>w[0]).slice(0,2).join('').toUpperCase();
const P=[['Byungju Kim','김병주','Ph.D., Physics','bjkim87@postech.ac.kr','Ph.D. Physics, POSTECH (2022) · B.S. Physics, Hanyang Univ. (2013)']];
const S=[
 ['Minsun Choi','최민선','Ph.D. candidate, Physics','mschoi1212@postech.ac.kr','B.S. Physics, Chungnam National Univ. (2024)'],
 ['Songmi Oh','오송미','Ph.D. candidate, IBIO','songmioh@postech.ac.kr','B.S. Biotechnology, Yonsei Univ. (2020)'],
 ['Philsang Yu','유필상','Ph.D. candidate, Physics','fx0908@postech.ac.kr','B.S. Physics, HUFS (2017)'],
 ['ManSeong Yoon','윤만성','Ph.D. candidate, Physics','yms5491@postech.ac.kr','B.S. Physics & Life Science, GIST (2025)'],
 ['Keunsang Yang','양근상','Ph.D. candidate, IBIO','yangks95@postech.ac.kr','B.S. Biotechnology, Konkuk Univ. (2018)'],
 ['GyeongHo Kang','강경호','Ph.D. candidate, Physics','ghkang@postech.ac.kr','B.S. Physics, GIST (2021)']];
const PHOTO={'Byungju Kim':'bjkim','Minsun Choi':'mschoi','Songmi Oh':'smoh','Philsang Yu':'psyu','ManSeong Yoon':'msyoon','Keunsang Yang':'ksyang','GyeongHo Kang':'ghkang'};
const card=p=>`<div class="person">${PHOTO[p[0]]?`<div class="avatar photo"><img src="img/people/${PHOTO[p[0]]}.webp" alt="${p[0]}" width="320" height="320" loading="lazy"></div>`:`<div class="avatar" aria-hidden="true">${init(p[0])}</div>`}<div><b>${p[0]} ${p[1]?`<span class="ko">(${p[1]})</span>`:''}</b><div class="deg">${p[2]}</div><div class="deg">${p[4]}</div><div class="mail">${p[3]}</div></div></div>`;
document.getElementById('postdocs').innerHTML=P.map(card).join('');
// graduate students in order of the year they received their B.S. (earliest first)
const bsYear=p=>+((p[4].match(/\((\d{4})\)/)||[])[1]||9999);
S.sort((a,b)=>bsYear(a)-bsYear(b));
document.getElementById('students').innerHTML=S.map(card).join('');
// [English name, Korean name, degree, year the degree was awarded, current position, homepage (optional, verified)]
const A=[['Jinseob Lee','이진섭','Ph.D.','2026',''],['Inho Yang','양인호','Ph.D.','2026',''],['Jincheol Seol','설진철','Ph.D.','2026',''],['Gayun Bu','부가연','Ph.D.','2024','Samsung Electronics'],['Jaeho Oh','오재호','MS','2024',''],['Byungju Kim','김병주','Ph.D.','2022','Postdoctoral Researcher, POSTECH'],['Daehyung Kim','김대형','Ph.D.','2019','Research Fellow, Boston Children\'s Hospital','https://labs.childrenshospital.org/kwonmoo-lee/people/daehyung-kim-phd'],['Ryanggeun Lee','이량근','Ph.D.','2019','Samsung Electronics'],['Minhyeok Chang','장민혁','Ph.D.','2018','Assistant Professor, KAIST','https://www.thechanglab.org/'],['Yongmoon Jeon','전용문','Ph.D.','2016','Samsung Display'],['Jong-Hyun Park','박종현','Ph.D.','2015','Lablup'],['Seung Hwan Lee','이승환','Ph.D.','2014','Associate Professor, Chung-Ang University','https://bio.cau.ac.kr/sub2_view.php?idx=100&tp='],['Sung-Sil Lee','이성실','Ph.D.','2013','Samsung Electronics'],['Won-Ki Cho','조원기','Ph.D.','2013','Associate Professor, KAIST','https://www.wonkicholab.com/'],['Soo-Ran Kim','김수란','MS','2011','Associate Professor, Kyungpook National University','https://sites.google.com/view/skimlab/home'],['Dong-Wook Kim','김동욱','MS','2009','Samsung Display'],['Daekil In','인대길','MS','2009','LG Electronics'],['Cherlhyun Jeong','','Postdoctoral Researcher','','KIST']];
// newest degree first; within a year, by last name, then first name
const surname=n=>n.trim().split(/\s+/).pop();
A.sort((x,y)=>y[3]-x[3]||surname(x[0]).localeCompare(surname(y[0]))||x[0].localeCompare(y[0]));
document.getElementById('alumni').innerHTML=A.map(([en,ko,deg,yr,now,url])=>`<div class="alum"><b>${url?`<a href="${url}">${en}</a>`:en}${ko?` <span class="ko">(${ko})</span>`:''}${url?` <a class="home" href="${url}" aria-label="${en} homepage">↗</a>`:''}</b><span class="deg">${[deg,yr].filter(Boolean).join(' ')}</span>${now?`<span class="now">${esc(now)}</span>`:''}</div>`).join('');

/* ---------- publications ---------- */
// [type, year, title, authors, venue, url, note]
const D=[
['paper',2026,'smDeepFLUOR: single-molecule deep learning fluorescence classification','Jinseob Lee, Byungju Kim, Gayun Bu, M. Tehseen, V.-S. Raducanu, S.M. Hamdan, Jong-Bong Lee*','Nature Communications 17, 7810','https://www.nature.com/articles/s41467-026-74716-3'],
['paper',2026,'Single-trajectory Bayesian modeling reveals multi-state diffusion of the MSH sliding clamp','S. Park#, Inho Yang#, Jinseob Lee, S. Kim, J. Martin-Lopez, R. Fishel, Jong-Bong Lee*, J.-H. Jeon*','npj Systems Biology and Applications 12, 98','https://www.nature.com/articles/s41540-026-00715-1'],
['paper',2026,'Tau condensation on DNA mediates microtubule attachment suggesting a mitotic role for centromere-localized tau','C. Park#, J. Jung#, Y. Hong#, H. Yoo, Keunsang Yang, J. Shin, M. Kim, C. Lim, A. Jeong, S. Hong, J.Y. Baek, S.-H. Rah, C. Lee-Eom, M. Seo, Y. Kim, J.-H. Jeon, Jong-Bong Lee, D.S. Hwang*, M.J. Shon*','Nature Communications 17, 1128','https://www.nature.com/articles/s41467-025-67888-x'],
['paper',2025,'Inferring spatially varying bending stiffness of biopolymers with deep-learning approach','C. Hong, C. Lim, Jaeho Oh, Jong-Bong Lee, E.E. Jeong, J.-Y. Yoo, W.K. Kim, J.-H. Jeon*','Machine Learning: Science and Technology 6, 045044','https://iopscience.iop.org/article/10.1088/2632-2153/ae1a35'],
['paper',2025,'Structural and dynamic basis of Ssp4-mediated DNA protection in foodborne bacterial spores','M. Seo, B. Kim, H. Shin, J. Kim, Jong-Bong Lee, Y.H. Ko, J.H. Kim, Y. Kim*','Scientific Reports 15, 37299','https://www.nature.com/articles/s41598-025-21159-3'],
['paper',2025,'DNA Hanger: surface-minimized single-molecule immunoassay platform','Jincheol Seol#, Byungju Kim#, E.-S. Yu, C. Jeong*, Jong-Bong Lee*','Small 21, 2409933','https://onlinelibrary.wiley.com/doi/10.1002/smll.202409933'],
['paper',2025,'A novel ADP-directed chaperone function facilitates the ATP-driven motor activity of SARS-CoV helicase','J. Yu, H. Im, H. Cho, Yongmoon Jeon, Jong-Bong Lee, G. Lee*','Nucleic Acids Research 53, gkaf034','https://academic.oup.com/nar/article/53/3/gkaf034/7986036'],
['paper',2025,'Attacking biological problems through single-molecule approaches','Jong-Bong Lee','BMB Reports 58(1), 1','https://www.bmbreports.org/journal/view.html?uid=2089&vmd=Full','Editorial'],
['paper',2025,'Cryogenic single-molecule fluorescence imaging','Phil Sang Yu#, C.U. Kim*, Jong-Bong Lee*','BMB Reports 58(1), 2–7','https://www.bmbreports.org/journal/view.html?uid=2090&vmd=Full'],
['paper',2023,'Single-molecule visualization of mRNA circularization during translation','Byungju Kim, Jincheol Seol, Y.K. Kim, Jong-Bong Lee*','Experimental & Molecular Medicine 55, 283–289','https://www.nature.com/articles/s12276-023-00933-1'],
['paper',2022,'Exploiting the distinctive properties of the bacterial and human MutS homolog sliding clamps on mismatched DNA','B.M. Britton, J.A. London, J. Martin-Lopez, N.D. Jones, J. Liu, Jong-Bong Lee, R. Fishel*','Journal of Biological Chemistry 298, 102505','https://www.sciencedirect.com/science/article/pii/S0021925822009486'],
['paper',2022,'Single polysome analysis of mRNP','Byungju Kim#, Y. Park#, H.J. Hwang, J. Chang, Y.K. Kim*, Jong-Bong Lee*','Biochemical and Biophysical Research Communications 618, 73–78','https://doi.org/10.1016/j.bbrc.2022.06.017'],
['paper',2022,'Formation of cellular close-ended tunneling nanotubes through mechanical deformation','Minhyeok Chang#, O-c Lee#, Gayun Bu#, Jaeho Oh, N.-O. Yunn, S.H. Ryu, H.-B. Kwon, A.B. Kolomeisky, S.-H. Shim, J. Doh, J.-H. Jeon*, Jong-Bong Lee*','Science Advances 8, eabj3995','https://www.science.org/doi/10.1126/sciadv.abj3995'],
['paper',2022,'The flip-flop configuration of the PABP-dimer leads to switching of the translation function','S. Gu, H.-M. Jeon, S.W. Nam, K.Y. Hong, M.S. Rahman, Jong-Bong Lee, Y. Kim, S.K. Jang*','Nucleic Acids Research 50, 306','https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab1205/6460780'],
['paper',2021,'TRIM28 functions as a negative regulator of aggresome formation','J. Chang, H.J. Hwang, Byungju Kim, Y.-G. Choi, J. Park, Y. Park, B.S. Lee, H. Park, M.J. Yoon, J.-S. Woo, C. Kim, M.-S. Park, Jong-Bong Lee, Y.K. Kim*','Autophagy','https://www.tandfonline.com/eprint/XEPCZK9RKUDWINF6S3IU/full?target=10.1080%2F15548627.2021.1909835'],
['paper',2021,'Linker domain function predicts pathogenic MLH1 missense variants','J.A. London#, J. Martin-Lopez#, Inho Yang#, J. Liu, Jong-Bong Lee*, R. Fishel*','Proceedings of the National Academy of Sciences of the USA 118(9), e2019215118','https://www.pnas.org/content/118/9/e2019215118'],
['paper',2020,'Nonsense-mediated mRNA decay factor UPF1 promotes aggresome formation','Y. Park#, J. Park#, H.J. Hwang#, Byungju Kim#, K. Jeong, J. Chang, Jong-Bong Lee*, Y.K. Kim*','Nature Communications 11, 3106','https://www.nature.com/articles/s41467-020-16939-6'],
['paper',2019,'MutL sliding clamps coordinate exonuclease-independent Escherichia coli mismatch repair','J. Liu#, Ryanggeun Lee#, B.M. Britton#, J.A. London, Keunsang Yang, J. Hanne, Jong-Bong Lee*, R. Fishel*','Nature Communications 10, 5294','https://www.nature.com/articles/s41467-019-13191-5'],
['paper',2019,'Multiplexed single-molecule flow-stretching bead assay for DNA enzymology','Ryanggeun Lee#, Keunsang Yang#, Jong-Bong Lee*','BMB Reports 52, 589','https://www.bmbreports.org/journal/view.html?uid=1416&vmd=Full&'],
['paper',2019,'Single-molecule DNA flow-stretching force spectroscopy with higher resolution using dark-field microscopy','Ryanggeun Lee, Yongmoon Jeon, Jong-Bong Lee*','Journal of the Korean Physical Society 75, 350','https://link.springer.com/article/10.3938/jkps.75.350'],
['paper',2019,'DNA skybridge: 3D structure producing a light sheet for high-throughput single-molecule imaging','Daehyung Kim, F. Rashid, Yeonmo Cho, M.S. Zaher, I.H. Cho, S.M. Hamdan, C. Jeong*, Jong-Bong Lee*','Nucleic Acids Research 47, e107','https://doi.org/10.1093/nar/gkz625'],
['paper',2018,'MutS homolog sliding clamps shield the DNA from binding proteins','J. Hanne, B.M. Britton, Jonghyun Park, J. Liu, J. Martin-Lopez, N. Jones, M. Schoffner, P. Klajner, R. Bundschuh*, Jong-Bong Lee*, R. Fishel*','Journal of Biological Chemistry 293, 14285','https://www.jbc.org/content/293/37/14285'],
['paper',2018,'Stochastic processes and component plasticity governing DNA mismatch repair','J. Liu, Jong-Bong Lee*, R. Fishel*','Journal of Molecular Biology 430, 4456','https://www.sciencedirect.com/science/article/pii/S0022283618305345'],
['paper',2018,'Coordinating multi-protein mismatch repair by managing diffusion mechanics on the DNA','Daehyung Kim, R. Fishel*, Jong-Bong Lee*','Journal of Molecular Biology 430, 4469','https://www.sciencedirect.com/science/article/pii/S0022283618304856'],
['paper',2018,'Dynamics of proofreading by the E. coli Pol III replicase','Jonghyun Park#, S. Jergic#, Yongmoon Jeon, Won-Ki Cho, Ryanggeun Lee, N.E. Dixon*, Jong-Bong Lee*','Cell Chemical Biology 25, 57','https://www.cell.com/cell-chemical-biology/fulltext/S2451-9456%2817%2930352-5'],
['paper',2017,'Extended depth of field for single biomolecule optical imaging-force spectroscopy','Minhyeok Chang, Jungsic Oh, Y. Kim, S. Hohng, Jong-Bong Lee*','Optics Express 25, 32189','https://www.osapublishing.org/oe/abstract.cfm?uri=oe-25-25-32189','OSA Spotlight · Nature Methods "Methods in Brief"'],
['paper',2017,'Light-induced fluorescence modulation of quantum dot–crystal violet conjugates: stochastic off–on–off cycles for multicolor patterning and super-resolution','S. Jung, J. Park, J. Bang, J.-Y. Kim, C. Kim, Yongmoon Jeon, Seung Hwan Lee, H. Jin, S. Choi, B. Kim, W.J. Lee, C.-G. Pack, Jong-Bong Lee, N.K. Lee, S. Kim*','Journal of the American Chemical Society 139, 7603','https://pubs.acs.org/doi/abs/10.1021/jacs.7b02530'],
['paper',2017,'The bent conformation of poly(A)-binding protein induced by RNA-binding is required for its translational activation function','K.Y. Hong, Seung Hwan Lee, S. Gu, E. Kim, S. An, J. Kwon, Jong-Bong Lee*, S.K. Jang*','RNA Biology 14, 370','https://www.tandfonline.com/doi/full/10.1080/15476286.2017.1280224'],
['paper',2016,'Cascading MutS and MutL sliding clamps control DNA diffusion to activate mismatch repair','J. Liu#, J. Hanne#, B.M. Britton, J. Bennett, Daehyung Kim, Jong-Bong Lee*, R. Fishel*','Nature 539, 583','https://www.nature.com/articles/nature20562','News & Views'],
['paper',2016,'Retroviral intasomes search for a target DNA by 1D diffusion which rarely results in integration','N.D. Jones, M.A. Lopez, J. Hanne, M.B. Peake, Jong-Bong Lee, R. Fishel*, K.E. Yoder*','Nature Communications 7, 11409','https://www.nature.com/ncomms/2016/160425/ncomms11409/full/ncomms11409.html'],
['paper',2016,'Pairwise detection of site-specific receptor phosphorylations using single-molecule blotting','K.L. Kim, Daehyung Kim, Seongsil Lee, S.-J. Kim, J.E. Noh, J.-H. Kim, Y.C. Chae, Jong-Bong Lee, S.H. Ryu*','Nature Communications 7, 11107','https://www.nature.com/articles/ncomms11107'],
['paper',2016,'Dynamic control of strand excision during human DNA mismatch repair','Yongmoon Jeon, Daehyung Kim, J.V. Martín-López, Ryanggeun Lee, Jungsic Oh, J. Hanne, R. Fishel*, Jong-Bong Lee*','Proceedings of the National Academy of Sciences of the USA 113, 3281','https://www.pnas.org/doi/abs/10.1073/pnas.1523748113','Commentary'],
['paper',2016,'Mismatch repair','R. Fishel, Jong-Bong Lee','DNA Replication, Recombination and Repair (Hanaoka & Sugasawa, eds.), Springer','https://link.springer.com/chapter/10.1007/978-4-431-55873-6_12','Book chapter'],
['paper',2015,'An efficient site-specific method for irreversible covalent labeling of proteins with a fluorophore','J. Liu, J. Hanne, B.M. Britton, M. Shoffner, A.E. Albers, J. Bennett, R. Zatezalo, R. Barfield, D. Rabuka, Jong-Bong Lee, R. Fishel*','Scientific Reports 5, 16883','https://www.nature.com/articles/srep16883'],
['paper',2015,'Widespread nuclease contamination in commonly used oxygen-scavenging systems','G. Senavirathne#, J. Liu#, M.A. Lopez Jr, J. Hanne, J. Martin-Lopez, Jong-Bong Lee, K.E. Yoder','Nature Methods 12, 901','https://www.nature.com/articles/nmeth.3588'],
['paper',2015,'Analysis of interactions between the epidermal growth factor receptor and soluble ligands on the basis of single-molecule diffusivity in the membrane of living cells','D.-H. Kim, K. Zhou, D.-K. Kim, S. Park, J. Noh, Y. Kwon, D. Kim, N.W. Song, Jong-Bong Lee, P.-G. Suh, N.K. Lee*, S.H. Ryu*','Angewandte Chemie International Edition 54, 7028','https://onlinelibrary.wiley.com/doi/full/10.1002/anie.201500871'],
['paper',2014,'Loading dynamics of a sliding DNA clamp','Won-Ki Cho, S. Jergic, Daehyung Kim, N.E. Dixon, Jong-Bong Lee*','Angewandte Chemie International Edition 53, 6768','https://onlinelibrary.wiley.com/doi/full/10.1002/anie.201403063'],
['paper',2014,'Crystal structures of the structure-selective nuclease Mus81-Eme1 bound to flap DNA substrates','G.H. Kwon, A. Jo, K. Baek, K.S. Jin, Y. Fu, Jong-Bong Lee, Y. Kim, Y. Cho*','The EMBO Journal 33, 1061','https://www.embopress.org/doi/full/10.1002/embj.201487820'],
['paper',2014,'Aptamer-based single-molecule imaging of insulin receptors in living cells','Minhyeok Chang#, Mijin Kwon#, Sooran Kim#, N.-O. Yunn#, Daehyung Kim, S.H. Ryu*, Jong-Bong Lee*','Journal of Biomedical Optics 19, 051204','https://www.spiedigitallibrary.org/journals/journal-of-biomedical-optics/volume-19/issue-05/051204/Aptamer-based-single-molecule-imaging-of-insulin-receptors-in-living/10.1117/1.JBO.19.5.051204.full'],
['paper',2014,'Single-molecule views of MutS on mismatched DNA','Jong-Bong Lee*, Won-Ki Cho, Jonghyun Park, Yongmoon Jeon, Seung Hwan Lee, R. Fishel*','DNA Repair 20, 82','https://www.sciencedirect.com/science/article/abs/pii/S1568786414000524','Invited review'],
['paper',2014,'Poly(A) RNA and Paip2 act as allosteric regulators of poly(A)-binding protein','Seung Hwan Lee, Jungsic Oh, Jonghyun Park, K.Y. Paek, S. Rho, S.K. Jang*, Jong-Bong Lee*','Nucleic Acids Research 42, 2697','https://academic.oup.com/nar/article/42/4/2697/2436199'],
['paper',2013,'Single-molecule FRET studies on DNA mismatch repair','J. Hanne, J. Liu, Jong-Bong Lee*, R. Fishel*','International Journal of Biophysics 3, 18','https://article.sapub.org/10.5923.s.biophysics.201311.03.html','Invited review'],
['paper',2013,'Mapping of surface-immobilized DNA with force-based atomic force microscopy','Y. Lee, S.H. Kwon, Y. Kim, Jong-Bong Lee, J.W. Park*','Analytical Chemistry 85, 4045','https://pubs.acs.org/doi/10.1021/ac3037848'],
['paper',2012,'Polarization-controlled photoswitching resolves dipole directions with subwavelength resolution','Seongsil Lee, Jungsic Oh, D. Kim, S. Kim, Jong-Bong Lee*, H.G. Nam*','Physical Review Letters 109, 248101','https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.248101'],
['paper',2012,'Assessment of radiation damage in single-shot coherent diffraction of DNA molecules by an extreme-ultraviolet free-electron laser','J. Park, D. Nam, Y. Kohmura, M. Nagasono, Yongmoon Jeon, Jong-Bong Lee, T. Ishikawa, C. Song*','Physical Review E 86, 042901','https://journals.aps.org/pre/abstract/10.1103/PhysRevE.86.042901'],
['paper',2012,'ATP alters the diffusion mechanics of MutS on mismatched DNA','Won-Ki Cho#, Cherlhyun Jeong#, Daehyung Kim, Minhyeok Chang, K.-M. Song, J. Hanne, C. Ban, R. Fishel*','Structure 20, 1264','https://www.cell.com/structure/fulltext/S0969-2126%2812%2900176-1'],
['paper',2011,'MutS switches between two fundamentally distinct clamps during mismatch repair','Cherlhyun Jeong#, Won-Ki Cho#, K.-M. Song, C. Cook, T.-Y. Yoon, C. Ban*, R. Fishel*, Jong-Bong Lee*','Nature Structural & Molecular Biology 18, 379','https://www.nature.com/articles/nsmb.2009'],
['paper',2011,'Unique temperature dependence and blinking behavior of CdTe/CdSe (core/shell) type-II quantum dots','B. Chon, Cherlhyun Jeong, J. Bang, J.H. Choi, Jong-Bong Lee, J. Park, T. Joo, S. Kim*','The Journal of Physical Chemistry C 115, 436','https://pubs.acs.org/doi/full/10.1021/jp109229u'],
['paper',2010,'Single-molecule analysis reveals the kinetics and physiological relevance of MutL-ssDNA binding','Daekil In, Jonghyun Park, R. Fishel*, C. Ban*, Jong-Bong Lee*','PLoS ONE 5(11), e15496','https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0015496'],
['paper',2007,'Dynamic DNA helicase–DNA polymerase interactions assure processive replication fork movement','A.M. van Oijen, D.E. Johnson, Jong-Bong Lee, N.A. Tanner, S.M. Hamdan, S. Tabor, U. Qimron','Molecular Cell 27, 539','https://www.sciencedirect.com/science/article/pii/S1097276507004091'],
['paper',2007,'Dynamics of the molecular orientation field coupled to ions in two-dimensional ferroelectric liquid crystals','R.A. Pelcovits, R.B. Meyer, Jong-Bong Lee*','Physical Review E 76, 021704','https://journals.aps.org/pre/abstract/10.1103/PhysRevE.76.021704'],
['paper',2007,'Role of electrostatics in the texture of islands in free-standing ferroelectric liquid crystal films','Jong-Bong Lee, R.A. Pelcovits, R.B. Meyer','Physical Review E 75, 051701','https://journals.aps.org/pre/abstract/10.1103/PhysRevE.75.051701'],
['paper',2006,'DNA primase acts as a molecular brake in DNA replication','Jong-Bong Lee, R.K. Hite, S.M. Hamdan, X.S. Xie, C.C. Richardson, A.M. van Oijen','Nature 439, 621','https://www.nature.com/articles/nature04317'],
['paper',2006,'Textural transformations in islands on free-standing smectic-C* liquid crystal films','Jong-Bong Lee, D. Konovalov, R.B. Meyer','Physical Review E 73, 051705','https://journals.aps.org/pre/abstract/10.1103/PhysRevE.73.051705'],
['paper',2001,'Equilibrium size and textures of islands in free-standing smectic C* films','R.B. Meyer, D. Konovalov, I. Kraus, Jong-Bong Lee','Molecular Crystals and Liquid Crystals 364, 123','https://www.tandfonline.com/doi/abs/10.1080/10587250108024982'],
['paper',1997,'Effect of the external twist on molecular orientation in antiferroelectric liquid crystals','J.-H. Lee, J.-H. Kim, Jong-Bong Lee, S.-D. Lee','Molecular Crystals and Liquid Crystals 301, 209','https://www.tandfonline.com/doi/abs/10.1080/10587259708041769'],
['pre',2025,'PAF15 stabilizes PCNA on DNA and slows down its sliding dynamics','Daehyung Kim, A. Gonzalez-Magaña, A.L. Chanu, Gayun Bu, M. Tehseen, C. Jeong, J.-H. Jeon, S.M. Hamdan, A. De Biasio*, F.J. Blanco*, Jong-Bong Lee*','bioRxiv 2025.11.15.688598','https://www.biorxiv.org/content/10.1101/2025.11.15.688598v1'],
['pre',2025,'smDeepFLUOR: single-molecule deep learning fluorescence classification','Jinseob Lee, Byungju Kim, Gayun Bu, M. Tehseen, V.-S. Raducanu, S.M. Hamdan, Jong-Bong Lee','bioRxiv 2025.07.09.663809','https://www.biorxiv.org/content/10.1101/2025.07.09.663809v1'],
['pre',2025,'Single-trajectory Bayesian modeling reveals multi-state diffusion of the MSH sliding clamp','S. Park, Inho Yang, Jinseob Lee, S. Kim, J. Martín-López, R. Fishel, Jong-Bong Lee, J.-H. Jeon','arXiv 2506.21943','https://arxiv.org/abs/2506.21943'],
['pre',2025,'Inferring spatially varying bending stiffness of biopolymers with deep-learning approach','C. Hong, C. Lim, Jaeho Oh, Jong-Bong Lee, E.E. Jeong, J.-Y. Yoo, W.K. Kim, J.-H. Jeon','bioRxiv 2025.06.12.659266','https://www.biorxiv.org/content/10.1101/2025.06.12.659266v1.abstract'],
['pre',2024,'DNA Hanger: surface free single-molecule blotting platform','Jincheol Seol, Byungju Kim, E.-S. Yu, C. Jeong, Jong-Bong Lee','bioRxiv 2024.10.23.619861','https://www.biorxiv.org/content/10.1101/2024.10.23.619861v1.abstract'],
['pre',2024,'Innovative strategies for enhancing signal-to-noise ratio in single-molecule imaging: the influence of molecular motion','Gayun Bu, Gyeongho Kang, S. Park, J.-H. Jeon, Jong-Bong Lee','bioRxiv 2024.10.13.618074','https://www.biorxiv.org/content/10.1101/2024.10.13.618074v1.abstract'],
['pre',2021,'The formation of a stable sliding clamp discriminates MSH2-MSH3 and MSH2-MSH6 mismatch interaction','B.M. Britton, J.A. London, J. Martin-Lopez, N.D. Jones, J. Liu, Jong-Bong Lee, R. Fishel','bioRxiv 2021.10.21.465318','https://www.biorxiv.org/content/10.1101/2021.10.21.465318v1.abstract'],
['pre',2018,'F-actin dynamics transform filopodial bridges into intercellular nanotubes capable of distant cell communication','Minhyeok Chang, Jaeho Oh, J. Doh, Jong-Bong Lee','bioRxiv 405340','https://www.biorxiv.org/content/10.1101/405340v1.abstract']
];
// [index in D, area class, area label, one-line significance]
const SEL=[
 [0,'tech','Techniques','A 3D convolutional network reads the emission patterns of single fluorophores to classify molecular states, including interactions of unlabeled proteins.'],
 [5,'tech','Techniques','A surface-minimized 3D scaffold suppresses nonspecific binding and reaches sub-picomolar detection in single-molecule immunoassays.'],
 [12,'cell','Cellular Dynamics','Proposes how close-ended tunneling nanotubes form: twisted double filopodial bridges turn into nanotubes under torsional and shearing stress.'],
 [15,'dna','DNA Mismatch Repair','Single-molecule analysis of the MLH1 linker domain distinguishes pathogenic missense variants found in Lynch syndrome.'],
 [16,'cell','Cellular Dynamics','The nonsense-mediated decay factor UPF1 is shown to promote aggresome formation.'],
 [17,'dna','DNA Mismatch Repair','MutL sliding clamps coordinate E. coli mismatch repair without the need for an exonuclease.'],
 [20,'tech','Techniques','DNA suspended on 3D quartz posts creates its own light sheet for high-throughput single-molecule imaging.'],
 [28,'dna','DNA Mismatch Repair','Shows MutS and MutL loading as cascading sliding clamps that control diffusion on DNA to switch on mismatch repair.'],
 [31,'dna','DNA Mismatch Repair','Real-time imaging shows how strand excision is dynamically controlled during human mismatch repair.'],
 [36,'dna','DNA Replication','Real-time view of how a sliding clamp is loaded onto DNA.'],
 [40,'cell','Translation','Poly(A) RNA and Paip2 act as allosteric regulators that change the conformation of poly(A)-binding protein.'],
 [43,'tech','Techniques','Polarization-controlled photoswitching resolves the dipole directions of single molecules below the diffraction limit.'],
 [46,'dna','DNA Mismatch Repair','smFRET reveals that MutS switches between two distinct clamp states during mismatch repair.'],
 [52,'dna','DNA Replication','DNA primase pauses the replication fork, acting as a molecular brake on leading-strand synthesis.']
];
const jname=v=>v.split(/ \d| \(/)[0];
document.getElementById('selected').innerHTML=SEL.map(([i,c,l,w])=>{const d=D[i];return `<li class="sel"><div class="jn">${esc(jname(d[4]))}<small>${d[1]} · ${esc(d[4].slice(jname(d[4]).length).trim())}</small></div><div class="body"><span class="area ${c}">${l}</span><a class="t" href="${d[5]}">${esc(d[2])}</a><p class="why">${esc(w)}</p><div class="au">${auth(d[3])}</div></div></li>`}).join('');
let filt='all',q='';
const list=document.getElementById('pub-list'),cnt=document.getElementById('pub-count');
function render(){
  const ql=q.trim().toLowerCase();
  const rows=D.filter(d=>(filt==='all'||d[0]===filt)&&(!ql||(d[2]+' '+d[3]+' '+d[4]).toLowerCase().includes(ql)));
  cnt.textContent=rows.length+' of '+D.length;
  if(!rows.length){list.innerHTML='<p class="empty">No publications match that search. Try an author surname or a journal name.</p>';return;}
  const by={};rows.forEach(d=>(by[d[1]]=by[d[1]]||[]).push(d));
  list.innerHTML=Object.keys(by).sort((a,b)=>b-a).map(y=>`<div class="year"><h3>${y}</h3><ol>${by[y].map(d=>`<li class="pub"><a class="t" href="${d[5]}">${esc(d[2])}</a>${d[0]==='pre'?'<span class="badge pre">preprint</span>':''}${d[6]?`<span class="badge note">${esc(d[6])}</span>`:''}<div class="au">${auth(d[3])}</div><div class="jr"><i>${esc(d[4])}</i> (${d[1]})</div></li>`).join('')}</ol></div>`).join('');
}
document.querySelectorAll('[data-f]').forEach(b=>b.addEventListener('click',()=>{filt=b.dataset.f;document.querySelectorAll('[data-f]').forEach(x=>x.setAttribute('aria-selected',x===b));render();}));
document.getElementById('pub-search').addEventListener('input',e=>{q=e.target.value;render();});
render();

/* ---------- copy ---------- */
document.querySelectorAll('.copy').forEach(b=>b.addEventListener('click',()=>{
  const done=()=>{b.textContent='Copied';setTimeout(()=>b.textContent='Copy',1500)};
  const sel=()=>{const s=b.previousElementSibling,r=document.createRange();r.selectNodeContents(s);const w=getSelection();w.removeAllRanges();w.addRange(r);b.textContent='Selected';setTimeout(()=>b.textContent='Copy',1500)};
  try{navigator.clipboard.writeText(b.dataset.v).then(done,sel)}catch(e){sel()}
}));

/* ---------- gallery ---------- */
// [title, year, posted date ('' for older albums), [[file, width, height], ...]] - converted from the old lab site
const GAL=[["Professor's birthday", 2026, "2026.08.14", [["48641_0", 1600, 1200]]], ["Inho and Jinseob's Graduation Celebration/Subin's Welcome Ceremony", 2026, "2026.07.29", [["48295_0", 1200, 1600], ["48295_1", 1200, 1600], ["48295_2", 1200, 1600], ["48295_3", 1200, 1600], ["48295_4", 1200, 1600], ["48295_5", 1200, 1600]]], ["Teacher's Day", 2026, "2026.05.11", [["46841_0", 1600, 1200], ["46841_1", 1600, 1200]]], ["Spring Group Photo", 2026, "2026.04.08", [["46164_0", 1274, 1600], ["46164_1", 1199, 1600]]], ["Jincheol's Graduation Ceremony", 2026, "2026.02.13", [["44498_0", 1600, 1200], ["44498_1", 1600, 1200], ["44498_2", 1600, 1078], ["44498_3", 1185, 1600]]], ["Jincheol's Graduation Celebration", 2026, "2026.01.13", [["44023_0", 1600, 1200], ["44023_1", 1600, 1200]]], ["Year End Party", 2025, "2025.12.30", [["43749_0", 1600, 1200], ["43749_1", 1600, 1200]]], ["LabMeeting/ Professor's birthday", 2025, "2025.08.13", [["28632_0", 1600, 1200], ["28632_1", 1600, 1200]]], ["Teacher's Day", 2025, "2025.06.17", [["23984_0", 1400, 1050], ["23983_0", 1553, 1165]]], ["Gayun Bu, Ph.D. 's Graduation Ceremony", 2025, "2025.06.17", [["23980_0", 1600, 1199], ["23979_0", 1400, 1050]]], ["Year End Party at Prof. Lee's Home", 2024, "", [["22767_0", 1600, 900]]], ["Yeongildae Park Picnic", 2023, "", [["23981_0", 1411, 1058]]], ["PhD thesis defense celebration for BJ Kim", 2022, "", [["22766_0", 1440, 1080]]], ["Teacher's Day", 2020, "", [["22765_0", 1600, 1200]]], ["Farewell Party for Dr.Kim", 2020, "", [["22764_0", 1280, 950]]], ["Birthday Party for Prof. Lee", 2020, "", [["22763_0", 1600, 1200]]], ["Birthday Party for Prof. Lee", 2019, "", [["22762_0", 1600, 1200]]], ["Doctoral Celebration for Dr. Lee & Dr. Kim", 2019, "", [["22761_0", 1024, 768]]], ["First Semester Closing Party", 2019, "", [["22760_0", 1024, 768]]], ["Teacher's Day", 2019, "", [["22759_0", 960, 721]]], ["Farewell party for Dr. Chang", 2019, "", [["22758_0", 1600, 1200]]], ["Christmas Party at Prof. Lee's Home", 2018, "", [["22757_0", 1024, 768]]], ["Second semester opening party", 2018, "", [["22756_0", 960, 720]]], ["Birthday Party for Keunsang", 2018, "", [["22755_0", 960, 720]]], ["Birthday Party for Prof. Lee", 2018, "", [["22754_0", 1024, 768]]], ["Workshop at Geoje Island", 2018, "", [["22753_0", 960, 720]]], ["First Semester Closing Party", 2018, "", [["22752_0", 960, 720]]], ["Teacher's Day", 2018, "", [["22751_0", 960, 720]]], ["Spring KPS with Prof. Jae-Hyung Jeon's group", 2018, "", [["22750_0", 960, 720]]], ["Lunch with Dr.Liu", 2018, "", [["22749_0", 960, 720]]], ["Year-end Party at Prof. Lee's home", 2017, "", [["22748_0", 960, 720]]], ["Lab 10th anniversary and homecoming day", 2017, "", [["22747_0", 960, 720]]], ["Birthday Party", 2017, "", [["22746_0", 960, 720]]], ["Teacher's Day", 2017, "", [["22745_0", 960, 720]]], ["Congratulations on Mueunjae-chaired Professor to Prof. Lee!", 2017, "", [["22744_0", 960, 640]]], ["Samsung HumanTech Awards, Minhyeok Chang(Gold Prize)", 2016, "", [["22743_0", 960, 640]]], ["Year-end Party", 2016, "", [["22742_0", 768, 576]]], ["August, farewell party for Dr. Jeon", 2016, "", [["22741_0", 768, 576], ["22740_0", 768, 576]]], ["Year-end party at Prof. Lee's House", 2015, "", [["22739_0", 960, 720]]], ["September, Birthday party", 2015, "", [["22738_0", 960, 720]]], ["Teacher's Day", 2015, "", [["22737_0", 960, 720]]], ["Teacher's Day", 2013, "", [["22736_0", 960, 720]]], ["Winter Workshop at the Vivaldi Park", 2013, "", [["22735_0", 960, 720]]], ["Year-end Party", 2012, "", [["22734_0", 960, 637]]], ["In the Snow", 2012, "", [["22733_0", 960, 720]]], ["New Year Party at Prof. Lee's House", 2012, "", [["22732_0", 960, 720]]], ["Summer outing. Sangju Seongjubong", 2011, "", [["22731_0", 1071, 1600]]], ["Winter outing. Phoenix Park", 2011, "", [["22730_0", 960, 640]]], ["Holiday party Prof. Lee's house", 2011, "", [["22729_0", 960, 720]]], ["Winter outing. Phoenix Park", 2010, "", [["22728_0", 960, 647]]], ["Teacher's Day", 2009, "", [["22727_0", 768, 514]]], ["Summer outing. Ha-Ok valley", 2009, "", [["22726_0", 1600, 522]]], ["Single-molecule group workshop", 2009, "", [["22723_0", 800, 246]]], ["Lab member", 2008, "", [["22725_0", 1024, 254]]], ["Lab member", 2007, "", [["22724_0", 1600, 412]]]];
const galList=document.getElementById('gal-list'),galYears=document.getElementById('gal-years'),galCount=document.getElementById('gal-count');
const years=[...new Set(GAL.map(a=>a[1]))];
let galYear='all',flat=[];
galYears.innerHTML=['all',...years].map(y=>`<button class="tab" type="button" data-y="${y}" aria-pressed="${y==='all'}">${y==='all'?'All':y}</button>`).join('');
function renderGal(){
  const sel=GAL.filter(a=>galYear==='all'||a[1]==galYear);
  flat=[];sel.forEach((a,ai)=>a[3].forEach((p,pi)=>flat.push({a,p,pi})));
  galCount.textContent=`${sel.length} albums \u00b7 ${flat.length} photos`;
  const by=[];sel.forEach(a=>{const g=by.find(x=>x.y===a[1]);g?g.a.push(a):by.push({y:a[1],a:[a]});});
  let k=0;
  galList.innerHTML=by.map(g=>`<h3 class="gal-year">${g.y}</h3><div class="gal-grid">${g.a.map(a=>{const start=k;k+=a[3].length;const c=a[3][0];
    return `<button class="album" type="button" data-i="${start}"><span class="cover"><img src="img/gallery/${c[0]}_t.webp" alt="" width="640" height="480" loading="lazy">${a[3].length>1?`<span class="count">${a[3].length} photos</span>`:''}</span><b>${esc(a[0])}</b><time>${a[2]||a[1]}</time></button>`}).join('')}</div>`).join('');
}
galYears.addEventListener('click',e=>{const b=e.target.closest('[data-y]');if(!b)return;galYear=b.dataset.y;galYears.querySelectorAll('[data-y]').forEach(x=>x.setAttribute('aria-pressed',x===b));renderGal();});
galList.addEventListener('click',e=>{const b=e.target.closest('.album');if(b)openLb(+b.dataset.i,b);});
renderGal();
const lb=document.getElementById('lb'),lbImg=document.getElementById('lb-img'),lbT=document.getElementById('lb-title'),lbS=document.getElementById('lb-sub');
const lbPrev=document.getElementById('lb-prev'),lbNext=document.getElementById('lb-next');
let lbI=0,lbFrom=null;
function showLb(){const it=flat[lbI];const [f,w,h]=it.p;lbImg.src=`img/gallery/${f}.webp`;lbImg.alt=it.a[0];
  lbT.textContent=it.a[0];lbS.textContent=`${it.a[2]||it.a[1]}${it.a[3].length>1?` \u00b7 ${it.pi+1} / ${it.a[3].length}`:''}`;
  lbPrev.disabled=lbI===0;lbNext.disabled=lbI===flat.length-1;
  [lbI-1,lbI+1].forEach(j=>{if(flat[j]){const im=new Image();im.src=`img/gallery/${flat[j].p[0]}.webp`;}});}
function openLb(i,from){lbI=i;lbFrom=from;lb.hidden=false;document.body.style.overflow='hidden';showLb();document.getElementById('lb-close').focus();}
function closeLb(){lb.hidden=true;document.body.style.overflow='';lbImg.removeAttribute('src');lbFrom&&lbFrom.focus();}
const step=d=>{const j=lbI+d;if(j>=0&&j<flat.length){lbI=j;showLb();}};
lbPrev.addEventListener('click',()=>step(-1));lbNext.addEventListener('click',()=>step(1));
document.getElementById('lb-close').addEventListener('click',closeLb);
lb.addEventListener('click',e=>{if(e.target===lb||e.target.classList.contains('lb-stage'))closeLb();});
addEventListener('keydown',e=>{if(lb.hidden)return;if(e.key==='Escape')closeLb();if(e.key==='ArrowLeft')step(-1);if(e.key==='ArrowRight')step(1);});
let tx=null;lb.addEventListener('touchstart',e=>{tx=e.touches[0].clientX},{passive:true});
lb.addEventListener('touchend',e=>{if(tx===null)return;const dx=e.changedTouches[0].clientX-tx;if(Math.abs(dx)>50)step(dx<0?1:-1);tx=null;});
addEventListener('hashchange',()=>{if(!lb.hidden)closeLb();});

/* ---------- home team preview ---------- */
const teamEl=document.getElementById('team-grid');
const tcard=(name,ko,role,photo,lead)=>`<li${lead?' class="lead"':''}>${photo?`<div class="avatar photo"><img src="img/people/${photo}.webp" alt="${name}" width="320" height="320" loading="lazy"></div>`:`<div class="avatar" aria-hidden="true">${init(name)}</div>`}<b>${name}</b><span class="ko">(${ko})</span><span class="role">${role}</span></li>`;
teamEl.innerHTML=tcard('Jong-Bong Lee','이종봉','Principal Investigator','jblee',true)
  +P.map(p=>tcard(p[0],p[1],'Postdoctoral Researcher',PHOTO[p[0]])).join('')
  +S.map(p=>tcard(p[0],p[1],'Ph.D. Candidate',PHOTO[p[0]])).join('');

/* ---------- views and navigation ---------- */
const VIEWS=['publications','members','gallery'];
const viewEls=[...document.querySelectorAll('.view')];
const navLinks=[...document.querySelectorAll('nav.main a')];
let current='';
function mark(key){navLinks.forEach(a=>{const on=a.dataset.view===key||a.dataset.sec===key;a.classList.toggle('on',on);on?a.setAttribute('aria-current','page'):a.removeAttribute('aria-current');});}
// keep section jumps clear of the sticky menu bar, whose height changes with screen width
const headerEl=document.querySelector('header.top');
const fitPad=()=>{document.documentElement.style.scrollPaddingTop=(headerEl.getBoundingClientRect().height+8)+'px';};
fitPad();addEventListener('resize',fitPad);if('ResizeObserver' in window)new ResizeObserver(fitPad).observe(headerEl);
const spyIds=['about','research','contact'];
let spyLock=0; // after a menu click, trust the click instead of the scroll position for a moment
function route(){
  let h=decodeURIComponent(location.hash.slice(1));
  if(h==='people'){history.replaceState(null,'','#members');h='members';} // old links
  const v=VIEWS.includes(h)?h:'home';
  if(v!==current){viewEls.forEach(el=>el.hidden=el.dataset.view!==v);current=v;}
  if(v!=='home'){window.scrollTo({top:0,behavior:'instant'});mark(v);document.title=v[0].toUpperCase()+v.slice(1)+' · JB Lab';return;}
  document.title='JB Lab \u00b7 Single-molecule Biophysics Lab, POSTECH';
  const target=h&&h!=='top'?document.getElementById(h):null;
  mark(spyIds.includes(h)?h:'');
  spyLock=performance.now()+900;
  if(target)requestAnimationFrame(()=>target.scrollIntoView({behavior:'instant',block:'start'}));
  else window.scrollTo({top:0,behavior:'instant'});
}
addEventListener('hashchange',route);
// clicking the link for the hash that is already in the address bar fires no hashchange, so route by hand
document.querySelectorAll('a[href^="#"]').forEach(a=>a.addEventListener('click',()=>{if(a.getAttribute('href')===location.hash)route();}));
route();

// on the home page, highlight the section that sits under the menu bar while scrolling
function spy(){
  if(current!=='home'||performance.now()<spyLock)return;
  const line=(document.querySelector('header.top')?.getBoundingClientRect().bottom||0)+48;
  let hit='';
  for(const id of ['about','research','highlights','team','contact']){
    const r=document.getElementById(id)?.getBoundingClientRect();
    if(r&&r.top<=line&&r.bottom>line){hit=id;break;}
  }
  // near the very bottom the short Contact section never reaches the line
  if(innerHeight+scrollY>=document.documentElement.scrollHeight-4)hit='contact';
  mark(spyIds.includes(hit)?hit:'');
}
let spyQueued=false;
addEventListener('scroll',()=>{if(!spyQueued){spyQueued=true;requestAnimationFrame(()=>{spyQueued=false;spy();});}},{passive:true});
})();

/* ---------- external links: new tab, no opener/referrer ---------- */
(function(){
  const fix=root=>root.querySelectorAll('a[href^="http"]').forEach(a=>{if(a.host!==location.host){a.target='_blank';a.rel='noopener noreferrer';}});
  fix(document);
  new MutationObserver(ms=>ms.forEach(m=>m.addedNodes.forEach(n=>{if(n.nodeType===1){if(n.matches('a[href^="http"]'))fix(n.parentNode);else fix(n);}}))).observe(document.body,{childList:true,subtree:true});
})();
