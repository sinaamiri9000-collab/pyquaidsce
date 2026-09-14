const fs=require('fs');
const vm=require('vm');
const html=fs.readFileSync('app/src/main/assets/index.html','utf8');
const m=html.match(/<script>([\s\S]*?)<\/script>/);
if(!m) throw new Error('embedded JS missing');

class FakeClassList{
  constructor(){this.s=new Set();}
  add(...x){x.forEach(v=>this.s.add(v));}
  remove(...x){x.forEach(v=>this.s.delete(v));}
  toggle(x,on){if(on===undefined){this.s.has(x)?this.s.delete(x):this.s.add(x);}else{on?this.s.add(x):this.s.delete(x);}}
  contains(x){return this.s.has(x);}
}
class FakeEl{
  constructor(){this.classList=new FakeClassList();this.style={};this.innerHTML='';this.textContent='';this.value='';this.disabled=false;this.dataset={};}
  addEventListener(){}
  closest(){return null;}
}
const els={};
function el(id){if(!els[id])els[id]=new FakeEl();return els[id];}
['content','loginForm','loginBtn','tokenInput','recentWrap','recentList','appView','loginView','profileRole','profileName','refreshBtn','toast','modal','modalTitle','modalBody'].forEach(el);
el('recentWrap').classList.add('hidden');
el('appView').classList.add('hidden');
const recent=[{token:'test-token-123456',role:'user',label:'حساب تست',lastUsed:Date.now()}];
const AndroidApi={
  request(action,payload){
    if(action==='resume')return JSON.stringify({ok:false});
    if(action==='recentTokens')return JSON.stringify({ok:true,data:recent});
    if(action==='logout')return JSON.stringify({ok:true,recentTokens:recent});
    if(action==='deactivate')return JSON.stringify({ok:true});
    return JSON.stringify({ok:true});
  },
  copy(){}, getDeviceId(){return 'test-device-id';}
};
const document={
  querySelector(sel){if(sel.startsWith('#'))return el(sel.slice(1));return new FakeEl();},
  querySelectorAll(){return[];},
  body:new FakeEl()
};
const ctx={console,document,AndroidApi,setTimeout:(fn)=>{fn();return 1;},clearTimeout(){},confirm:()=>true,Intl,Date,Math,JSON,String,Number,Array,Object,Error,RegExp,encodeURIComponent,decodeURIComponent};
vm.createContext(ctx);
vm.runInContext(m[1],ctx,{timeout:3000});

// Recent tokens must be visible immediately, including immediately after logout.
if(el('recentWrap').classList.contains('hidden')) throw new Error('recent tokens hidden on login screen');
vm.runInContext("state={role:'user',data:{},token:'x',activeSubscriptionId:0,activeAgentType:'',activeTarget:''}; logout();",ctx);
if(el('recentWrap').classList.contains('hidden')) throw new Error('recent tokens hidden after logout');
if(!el('recentList').innerHTML.includes('حساب تست')) throw new Error('recent token row not rendered after logout');

// Active subscription UI mirrors extension semantics: activation button is hidden, login/deactivate are visible.
let active=vm.runInContext("state.activeSubscriptionId=12; serviceHtml({type:'GPT',subscriptions:[{id:12,account:{name:'A'}}]})",ctx);
if(!active.includes('ورود به سرویس')||!active.includes('غیرفعال‌سازی')) throw new Error('active subscription controls missing');
if(active.includes('فعال‌سازی و ورود')) throw new Error('activation button visible for active subscription');
let inactive=vm.runInContext("state.activeSubscriptionId=0; serviceHtml({type:'GPT',subscriptions:[{id:12,account:{name:'A'}}]})",ctx);
if(!inactive.includes('فعال‌سازی و ورود')) throw new Error('inactive activation control missing');

// Jalali conversion regression.
let g=vm.runInContext('toGregorian(1405,6,23)',ctx);
let j=vm.runInContext(`toJalaali(${g.gy},${g.gm},${g.gd})`,ctx);
if(j.jy!==1405||j.jm!==6||j.jd!==23) throw new Error('Jalali round-trip failed');
console.log('v1.4 UI regression tests passed');
