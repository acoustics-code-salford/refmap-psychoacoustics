"use strict";(self.webpackChunk_lckr_jupyterlab_variableinspector=self.webpackChunk_lckr_jupyterlab_variableinspector||[]).push([[372],{3180:(r,t,e)=>{e.d(t,{$2:()=>k,DR:()=>v,KW:()=>x,T8:()=>w,hM:()=>h,hP:()=>f,iI:()=>p,lw:()=>d,rD:()=>y,rp:()=>u,v1:()=>N,wo:()=>c,zP:()=>m});var n=e(8094),i=e(9891),s=e(906),a=e(6874),b=e(1404),g=e(5400),o=e(5688);function u(r){return.2126*r.r+.7152*r.g+.0722*r.b}function h(r){function t(r){return r<=.03928?r/12.92:Math.pow((r+.055)/1.055,2.4)}return u(new b.h(t(r.r),t(r.g),t(r.b),1))}const l=(r,t)=>(r+.05)/(t+.05);function c(r,t){const e=h(r),n=h(t);return e>n?l(e,n):l(n,e)}function d(r){const t=Math.max(r.r,r.g,r.b),e=Math.min(r.r,r.g,r.b),i=t-e;let s=0;0!==i&&(s=t===r.r?(r.g-r.b)/i%6*60:t===r.g?60*((r.b-r.r)/i+2):60*((r.r-r.g)/i+4)),s<0&&(s+=360);const a=(t+e)/2;let b=0;return 0!==i&&(b=i/(1-Math.abs(2*a-1))),new n.H(s,b,a)}function f(r,t=1){const e=(1-Math.abs(2*r.l-1))*r.s,n=e*(1-Math.abs(r.h/60%2-1)),i=r.l-e/2;let s=0,a=0,g=0;return r.h<60?(s=e,a=n,g=0):r.h<120?(s=n,a=e,g=0):r.h<180?(s=0,a=e,g=n):r.h<240?(s=0,a=n,g=e):r.h<300?(s=n,a=0,g=e):r.h<360&&(s=e,a=0,g=n),new b.h(s+i,a+i,g+i,t)}function w(r){const t=Math.max(r.r,r.g,r.b),e=t-Math.min(r.r,r.g,r.b);let n=0;0!==e&&(n=t===r.r?(r.g-r.b)/e%6*60:t===r.g?60*((r.b-r.r)/e+2):60*((r.r-r.g)/e+4)),n<0&&(n+=360);let s=0;return 0!==t&&(s=e/t),new i.T(n,s,t)}function p(r,t=1){const e=r.s*r.v,n=e*(1-Math.abs(r.h/60%2-1)),i=r.v-e;let s=0,a=0,g=0;return r.h<60?(s=e,a=n,g=0):r.h<120?(s=n,a=e,g=0):r.h<180?(s=0,a=e,g=n):r.h<240?(s=0,a=n,g=e):r.h<300?(s=n,a=0,g=e):r.h<360&&(s=e,a=0,g=n),new b.h(s+i,a+i,g+i,t)}function m(r){function t(r){return r<=.04045?r/12.92:Math.pow((r+.055)/1.055,2.4)}const e=t(r.r),n=t(r.g),i=t(r.b),s=.4124564*e+.3575761*n+.1804375*i,a=.2126729*e+.7151522*n+.072175*i,b=.0193339*e+.119192*n+.9503041*i;return new g.x(s,a,b)}function y(r,t=1){function e(r){return r<=.0031308?12.92*r:1.055*Math.pow(r,1/2.4)-.055}const n=e(3.2404542*r.x-1.5371385*r.y-.4985314*r.z),i=e(-.969266*r.x+1.8760108*r.y+.041556*r.z),s=e(.0556434*r.x-.2040259*r.y+1.0572252*r.z);return new b.h(n,i,s,t)}function N(r){return function(r){function t(r){return r>s.R.epsilon?Math.pow(r,1/3):(s.R.kappa*r+16)/116}const e=t(r.x/g.x.whitePoint.x),n=t(r.y/g.x.whitePoint.y),i=116*n-16,a=500*(e-n),b=200*(n-t(r.z/g.x.whitePoint.z));return new s.R(i,a,b)}(m(r))}function v(r,t=1){return y(function(r){const t=(r.l+16)/116,e=t+r.a/500,n=t-r.b/200,i=Math.pow(e,3),a=Math.pow(t,3),b=Math.pow(n,3);let o=0;o=i>s.R.epsilon?i:(116*e-16)/s.R.kappa;let u=0;u=r.l>s.R.epsilon*s.R.kappa?a:r.l/s.R.kappa;let h=0;return h=b>s.R.epsilon?b:(116*n-16)/s.R.kappa,o=g.x.whitePoint.x*o,u=g.x.whitePoint.y*u,h=g.x.whitePoint.z*h,new g.x(o,u,h)}(r),t)}function k(r){return function(r){let t=0;(Math.abs(r.b)>.001||Math.abs(r.a)>.001)&&(t=(0,o.vi)(Math.atan2(r.b,r.a))),t<0&&(t+=360);const e=Math.sqrt(r.a*r.a+r.b*r.b);return new a.t(r.l,e,t)}(N(r))}function x(r,t=1){return v(function(r){let t=0,e=0;return 0!==r.h&&(t=Math.cos((0,o.Ht)(r.h))*r.c,e=Math.sin((0,o.Ht)(r.h))*r.c),new s.R(r.l,t,e)}(r),t)}},8094:(r,t,e)=>{e.d(t,{H:()=>i});var n=e(5688);class i{constructor(r,t,e){this.h=r,this.s=t,this.l=e}static fromObject(r){return!r||isNaN(r.h)||isNaN(r.s)||isNaN(r.l)?null:new i(r.h,r.s,r.l)}equalValue(r){return this.h===r.h&&this.s===r.s&&this.l===r.l}roundToPrecision(r){return new i((0,n.fZ)(this.h,r),(0,n.fZ)(this.s,r),(0,n.fZ)(this.l,r))}toObject(){return{h:this.h,s:this.s,l:this.l}}}},9891:(r,t,e)=>{e.d(t,{T:()=>i});var n=e(5688);class i{constructor(r,t,e){this.h=r,this.s=t,this.v=e}static fromObject(r){return!r||isNaN(r.h)||isNaN(r.s)||isNaN(r.v)?null:new i(r.h,r.s,r.v)}equalValue(r){return this.h===r.h&&this.s===r.s&&this.v===r.v}roundToPrecision(r){return new i((0,n.fZ)(this.h,r),(0,n.fZ)(this.s,r),(0,n.fZ)(this.v,r))}toObject(){return{h:this.h,s:this.s,v:this.v}}}},906:(r,t,e)=>{e.d(t,{R:()=>i});var n=e(5688);class i{constructor(r,t,e){this.l=r,this.a=t,this.b=e}static fromObject(r){return!r||isNaN(r.l)||isNaN(r.a)||isNaN(r.b)?null:new i(r.l,r.a,r.b)}equalValue(r){return this.l===r.l&&this.a===r.a&&this.b===r.b}roundToPrecision(r){return new i((0,n.fZ)(this.l,r),(0,n.fZ)(this.a,r),(0,n.fZ)(this.b,r))}toObject(){return{l:this.l,a:this.a,b:this.b}}}i.epsilon=216/24389,i.kappa=24389/27},6874:(r,t,e)=>{e.d(t,{t:()=>i});var n=e(5688);class i{constructor(r,t,e){this.l=r,this.c=t,this.h=e}static fromObject(r){return!r||isNaN(r.l)||isNaN(r.c)||isNaN(r.h)?null:new i(r.l,r.c,r.h)}equalValue(r){return this.l===r.l&&this.c===r.c&&this.h===r.h}roundToPrecision(r){return new i((0,n.fZ)(this.l,r),(0,n.fZ)(this.c,r),(0,n.fZ)(this.h,r))}toObject(){return{l:this.l,c:this.c,h:this.h}}}},1404:(r,t,e)=>{e.d(t,{h:()=>i});var n=e(5688);class i{constructor(r,t,e,n){this.r=r,this.g=t,this.b=e,this.a="number"!=typeof n||isNaN(n)?1:n}static fromObject(r){return!r||isNaN(r.r)||isNaN(r.g)||isNaN(r.b)?null:new i(r.r,r.g,r.b,r.a)}equalValue(r){return this.r===r.r&&this.g===r.g&&this.b===r.b&&this.a===r.a}toStringHexRGB(){return"#"+[this.r,this.g,this.b].map(this.formatHexValue).join("")}toStringHexRGBA(){return this.toStringHexRGB()+this.formatHexValue(this.a)}toStringHexARGB(){return"#"+[this.a,this.r,this.g,this.b].map(this.formatHexValue).join("")}toStringWebRGB(){return`rgb(${Math.round((0,n.cY)(this.r,0,255))},${Math.round((0,n.cY)(this.g,0,255))},${Math.round((0,n.cY)(this.b,0,255))})`}toStringWebRGBA(){return`rgba(${Math.round((0,n.cY)(this.r,0,255))},${Math.round((0,n.cY)(this.g,0,255))},${Math.round((0,n.cY)(this.b,0,255))},${(0,n.uZ)(this.a,0,1)})`}roundToPrecision(r){return new i((0,n.fZ)(this.r,r),(0,n.fZ)(this.g,r),(0,n.fZ)(this.b,r),(0,n.fZ)(this.a,r))}clamp(){return new i((0,n.uZ)(this.r,0,1),(0,n.uZ)(this.g,0,1),(0,n.uZ)(this.b,0,1),(0,n.uZ)(this.a,0,1))}toObject(){return{r:this.r,g:this.g,b:this.b,a:this.a}}formatHexValue(r){return(0,n.yi)((0,n.cY)(r,0,255))}}},5400:(r,t,e)=>{e.d(t,{x:()=>i});var n=e(5688);class i{constructor(r,t,e){this.x=r,this.y=t,this.z=e}static fromObject(r){return!r||isNaN(r.x)||isNaN(r.y)||isNaN(r.z)?null:new i(r.x,r.y,r.z)}equalValue(r){return this.x===r.x&&this.y===r.y&&this.z===r.z}roundToPrecision(r){return new i((0,n.fZ)(this.x,r),(0,n.fZ)(this.y,r),(0,n.fZ)(this.z,r))}toObject(){return{x:this.x,y:this.y,z:this.z}}}i.whitePoint=new i(.95047,1,1.08883)},5688:(r,t,e)=>{function n(r,t,e){return isNaN(r)||r<=t?t:r>=e?e:r}function i(r,t,e){return isNaN(r)||r<=t?0:r>=e?1:r/(e-t)}function s(r,t,e){return isNaN(r)?t:t+r*(e-t)}function a(r){return r*(Math.PI/180)}function b(r){return r*(180/Math.PI)}function g(r){const t=Math.round(n(r,0,255)).toString(16);return 1===t.length?"0"+t:t}function o(r,t,e){return isNaN(r)||r<=0?t:r>=1?e:t+r*(e-t)}function u(r,t,e){if(r<=0)return t%360;if(r>=1)return e%360;const n=(t-e+360)%360;return n<=(e-t+360)%360?(t-n*r+360)%360:(t+n*r+360)%360}function h(r,t){const e=Math.pow(10,t);return Math.round(r*e)/e}e.d(t,{AG:()=>u,Fv:()=>i,Ht:()=>a,cY:()=>s,fZ:()=>h,t7:()=>o,uZ:()=>n,vi:()=>b,yi:()=>g}),Math.PI},7958:(r,t,e)=>{e.d(t,{lu:()=>h,in:()=>u});var n=e(1404),i=e(5688);const s={aliceblue:{r:.941176,g:.972549,b:1},antiquewhite:{r:.980392,g:.921569,b:.843137},aqua:{r:0,g:1,b:1},aquamarine:{r:.498039,g:1,b:.831373},azure:{r:.941176,g:1,b:1},beige:{r:.960784,g:.960784,b:.862745},bisque:{r:1,g:.894118,b:.768627},black:{r:0,g:0,b:0},blanchedalmond:{r:1,g:.921569,b:.803922},blue:{r:0,g:0,b:1},blueviolet:{r:.541176,g:.168627,b:.886275},brown:{r:.647059,g:.164706,b:.164706},burlywood:{r:.870588,g:.721569,b:.529412},cadetblue:{r:.372549,g:.619608,b:.627451},chartreuse:{r:.498039,g:1,b:0},chocolate:{r:.823529,g:.411765,b:.117647},coral:{r:1,g:.498039,b:.313725},cornflowerblue:{r:.392157,g:.584314,b:.929412},cornsilk:{r:1,g:.972549,b:.862745},crimson:{r:.862745,g:.078431,b:.235294},cyan:{r:0,g:1,b:1},darkblue:{r:0,g:0,b:.545098},darkcyan:{r:0,g:.545098,b:.545098},darkgoldenrod:{r:.721569,g:.52549,b:.043137},darkgray:{r:.662745,g:.662745,b:.662745},darkgreen:{r:0,g:.392157,b:0},darkgrey:{r:.662745,g:.662745,b:.662745},darkkhaki:{r:.741176,g:.717647,b:.419608},darkmagenta:{r:.545098,g:0,b:.545098},darkolivegreen:{r:.333333,g:.419608,b:.184314},darkorange:{r:1,g:.54902,b:0},darkorchid:{r:.6,g:.196078,b:.8},darkred:{r:.545098,g:0,b:0},darksalmon:{r:.913725,g:.588235,b:.478431},darkseagreen:{r:.560784,g:.737255,b:.560784},darkslateblue:{r:.282353,g:.239216,b:.545098},darkslategray:{r:.184314,g:.309804,b:.309804},darkslategrey:{r:.184314,g:.309804,b:.309804},darkturquoise:{r:0,g:.807843,b:.819608},darkviolet:{r:.580392,g:0,b:.827451},deeppink:{r:1,g:.078431,b:.576471},deepskyblue:{r:0,g:.74902,b:1},dimgray:{r:.411765,g:.411765,b:.411765},dimgrey:{r:.411765,g:.411765,b:.411765},dodgerblue:{r:.117647,g:.564706,b:1},firebrick:{r:.698039,g:.133333,b:.133333},floralwhite:{r:1,g:.980392,b:.941176},forestgreen:{r:.133333,g:.545098,b:.133333},fuchsia:{r:1,g:0,b:1},gainsboro:{r:.862745,g:.862745,b:.862745},ghostwhite:{r:.972549,g:.972549,b:1},gold:{r:1,g:.843137,b:0},goldenrod:{r:.854902,g:.647059,b:.12549},gray:{r:.501961,g:.501961,b:.501961},green:{r:0,g:.501961,b:0},greenyellow:{r:.678431,g:1,b:.184314},grey:{r:.501961,g:.501961,b:.501961},honeydew:{r:.941176,g:1,b:.941176},hotpink:{r:1,g:.411765,b:.705882},indianred:{r:.803922,g:.360784,b:.360784},indigo:{r:.294118,g:0,b:.509804},ivory:{r:1,g:1,b:.941176},khaki:{r:.941176,g:.901961,b:.54902},lavender:{r:.901961,g:.901961,b:.980392},lavenderblush:{r:1,g:.941176,b:.960784},lawngreen:{r:.486275,g:.988235,b:0},lemonchiffon:{r:1,g:.980392,b:.803922},lightblue:{r:.678431,g:.847059,b:.901961},lightcoral:{r:.941176,g:.501961,b:.501961},lightcyan:{r:.878431,g:1,b:1},lightgoldenrodyellow:{r:.980392,g:.980392,b:.823529},lightgray:{r:.827451,g:.827451,b:.827451},lightgreen:{r:.564706,g:.933333,b:.564706},lightgrey:{r:.827451,g:.827451,b:.827451},lightpink:{r:1,g:.713725,b:.756863},lightsalmon:{r:1,g:.627451,b:.478431},lightseagreen:{r:.12549,g:.698039,b:.666667},lightskyblue:{r:.529412,g:.807843,b:.980392},lightslategray:{r:.466667,g:.533333,b:.6},lightslategrey:{r:.466667,g:.533333,b:.6},lightsteelblue:{r:.690196,g:.768627,b:.870588},lightyellow:{r:1,g:1,b:.878431},lime:{r:0,g:1,b:0},limegreen:{r:.196078,g:.803922,b:.196078},linen:{r:.980392,g:.941176,b:.901961},magenta:{r:1,g:0,b:1},maroon:{r:.501961,g:0,b:0},mediumaquamarine:{r:.4,g:.803922,b:.666667},mediumblue:{r:0,g:0,b:.803922},mediumorchid:{r:.729412,g:.333333,b:.827451},mediumpurple:{r:.576471,g:.439216,b:.858824},mediumseagreen:{r:.235294,g:.701961,b:.443137},mediumslateblue:{r:.482353,g:.407843,b:.933333},mediumspringgreen:{r:0,g:.980392,b:.603922},mediumturquoise:{r:.282353,g:.819608,b:.8},mediumvioletred:{r:.780392,g:.082353,b:.521569},midnightblue:{r:.098039,g:.098039,b:.439216},mintcream:{r:.960784,g:1,b:.980392},mistyrose:{r:1,g:.894118,b:.882353},moccasin:{r:1,g:.894118,b:.709804},navajowhite:{r:1,g:.870588,b:.678431},navy:{r:0,g:0,b:.501961},oldlace:{r:.992157,g:.960784,b:.901961},olive:{r:.501961,g:.501961,b:0},olivedrab:{r:.419608,g:.556863,b:.137255},orange:{r:1,g:.647059,b:0},orangered:{r:1,g:.270588,b:0},orchid:{r:.854902,g:.439216,b:.839216},palegoldenrod:{r:.933333,g:.909804,b:.666667},palegreen:{r:.596078,g:.984314,b:.596078},paleturquoise:{r:.686275,g:.933333,b:.933333},palevioletred:{r:.858824,g:.439216,b:.576471},papayawhip:{r:1,g:.937255,b:.835294},peachpuff:{r:1,g:.854902,b:.72549},peru:{r:.803922,g:.521569,b:.247059},pink:{r:1,g:.752941,b:.796078},plum:{r:.866667,g:.627451,b:.866667},powderblue:{r:.690196,g:.878431,b:.901961},purple:{r:.501961,g:0,b:.501961},red:{r:1,g:0,b:0},rosybrown:{r:.737255,g:.560784,b:.560784},royalblue:{r:.254902,g:.411765,b:.882353},saddlebrown:{r:.545098,g:.270588,b:.07451},salmon:{r:.980392,g:.501961,b:.447059},sandybrown:{r:.956863,g:.643137,b:.376471},seagreen:{r:.180392,g:.545098,b:.341176},seashell:{r:1,g:.960784,b:.933333},sienna:{r:.627451,g:.321569,b:.176471},silver:{r:.752941,g:.752941,b:.752941},skyblue:{r:.529412,g:.807843,b:.921569},slateblue:{r:.415686,g:.352941,b:.803922},slategray:{r:.439216,g:.501961,b:.564706},slategrey:{r:.439216,g:.501961,b:.564706},snow:{r:1,g:.980392,b:.980392},springgreen:{r:0,g:1,b:.498039},steelblue:{r:.27451,g:.509804,b:.705882},tan:{r:.823529,g:.705882,b:.54902},teal:{r:0,g:.501961,b:.501961},thistle:{r:.847059,g:.74902,b:.847059},tomato:{r:1,g:.388235,b:.278431},transparent:{r:0,g:0,b:0,a:0},turquoise:{r:.25098,g:.878431,b:.815686},violet:{r:.933333,g:.509804,b:.933333},wheat:{r:.960784,g:.870588,b:.701961},white:{r:1,g:1,b:1},whitesmoke:{r:.960784,g:.960784,b:.960784},yellow:{r:1,g:1,b:0},yellowgreen:{r:.603922,g:.803922,b:.196078}},a=/^rgb\(\s*((?:(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*,\s*){2}(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*)\)$/i,b=/^rgba\(\s*((?:(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*,\s*){3}(?:0|1|0?\.\d*)\s*)\)$/i,g=/^#((?:[0-9a-f]{6}|[0-9a-f]{3}))$/i,o=/^#((?:[0-9a-f]{8}|[0-9a-f]{4}))$/i;function u(r){const t=g.exec(r);if(null===t)return null;let e=t[1];if(3===e.length){const r=e.charAt(0),t=e.charAt(1),n=e.charAt(2);e=r.concat(r,t,t,n,n)}const s=parseInt(e,16);return isNaN(s)?null:new n.h((0,i.Fv)((16711680&s)>>>16,0,255),(0,i.Fv)((65280&s)>>>8,0,255),(0,i.Fv)(255&s,0,255),1)}function h(r){const t=r.toLowerCase();return function(r){return g.test(r)}(t)?u(t):function(r){return function(r){return o.test(r)}(r)}(t)?function(r){const t=o.exec(r);if(null===t)return null;let e=t[1];if(4===e.length){const r=e.charAt(0),t=e.charAt(1),n=e.charAt(2),i=e.charAt(3);e=r.concat(r,t,t,n,n,i,i)}const s=parseInt(e,16);return isNaN(s)?null:new n.h((0,i.Fv)((16711680&s)>>>16,0,255),(0,i.Fv)((65280&s)>>>8,0,255),(0,i.Fv)(255&s,0,255),(0,i.Fv)((4278190080&s)>>>24,0,255))}(t):function(r){return a.test(r)}(t)?function(r){const t=a.exec(r);if(null===t)return null;const e=t[1].split(",");return new n.h((0,i.Fv)(Number(e[0]),0,255),(0,i.Fv)(Number(e[1]),0,255),(0,i.Fv)(Number(e[2]),0,255),1)}(t):function(r){return b.test(r)}(t)?function(r){const t=b.exec(r);if(null===t)return null;const e=t[1].split(",");return 4===e.length?new n.h((0,i.Fv)(Number(e[0]),0,255),(0,i.Fv)(Number(e[1]),0,255),(0,i.Fv)(Number(e[2]),0,255),Number(e[3])):null}(t):function(r){return s.hasOwnProperty(r)}(t)?function(r){const t=s[r.toLowerCase()];return t?new n.h(t.r,t.g,t.b,t.hasOwnProperty("a")?t.a:void 0):null}(t):null}},655:(r,t,e)=>{function n(r,t,e,n){var i,s=arguments.length,a=s<3?t:null===n?n=Object.getOwnPropertyDescriptor(t,e):n;if("object"==typeof Reflect&&"function"==typeof Reflect.decorate)a=Reflect.decorate(r,t,e,n);else for(var b=r.length-1;b>=0;b--)(i=r[b])&&(a=(s<3?i(a):s>3?i(t,e,a):i(t,e))||a);return s>3&&a&&Object.defineProperty(t,e,a),a}e.d(t,{gn:()=>n})}}]);