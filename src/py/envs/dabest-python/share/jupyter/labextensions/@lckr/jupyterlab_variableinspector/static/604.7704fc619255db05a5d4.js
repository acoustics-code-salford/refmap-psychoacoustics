/*! For license information please see 604.7704fc619255db05a5d4.js.LICENSE.txt */
"use strict";(self.webpackChunk_lckr_jupyterlab_variableinspector=self.webpackChunk_lckr_jupyterlab_variableinspector||[]).push([[604],{604:(t,e,i)=>{i.r(e),i.d(e,{ARIAGlobalStatesAndProperties:()=>Rt,Accordion:()=>Ot,AccordionExpandMode:()=>Tt,AccordionItem:()=>ct,Anchor:()=>St,AnchoredRegion:()=>qt,Avatar:()=>Qt,Badge:()=>Jt,BaseProgress:()=>Hs,Breadcrumb:()=>se,BreadcrumbItem:()=>ee,Button:()=>pe,Calendar:()=>ge,CalendarTitleTemplate:()=>ke,Card:()=>Fe,CheckableFormAssociated:()=>de,Checkbox:()=>He,Combobox:()=>Qe,ComboboxAutocomplete:()=>Xe,ComponentPresentation:()=>nt,Container:()=>C,ContainerConfiguration:()=>v,ContainerImpl:()=>K,DI:()=>y,DataGrid:()=>Ie,DataGridCell:()=>we,DataGridCellTypes:()=>fe,DataGridRow:()=>$e,DataGridRowTypes:()=>ye,DateFormatter:()=>ve,DefaultComponentPresentation:()=>at,DefaultResolver:()=>m,DelegatesARIAButton:()=>me,DelegatesARIACombobox:()=>Ze,DelegatesARIALink:()=>At,DelegatesARIAListbox:()=>Ke,DelegatesARIAListboxOption:()=>je,DelegatesARIASearch:()=>Zs,DelegatesARIASelect:()=>io,DelegatesARIATextbox:()=>As,DelegatesARIAToolbar:()=>Mo,DesignSystem:()=>Ti,DesignToken:()=>xi,Dialog:()=>Ui,Disclosure:()=>Yi,Divider:()=>Ji,DividerRole:()=>Zi,ElementDisambiguation:()=>wi,FactoryImpl:()=>V,Flipper:()=>is,FlipperDirection:()=>ts,FlyoutPosBottom:()=>_t,FlyoutPosBottomFill:()=>Gt,FlyoutPosTallest:()=>Kt,FlyoutPosTallestFill:()=>Yt,FlyoutPosTop:()=>jt,FlyoutPosTopFill:()=>Wt,FormAssociated:()=>he,FoundationElement:()=>rt,FoundationElementRegistry:()=>ht,GenerateHeaderOptions:()=>be,HorizontalScroll:()=>_s,Listbox:()=>_e,ListboxElement:()=>os,ListboxOption:()=>Ue,MatchMediaBehavior:()=>_o,MatchMediaStyleSheetBehavior:()=>Ko,Menu:()=>Es,MenuItem:()=>Is,MenuItemRole:()=>xs,NumberField:()=>Ps,Picker:()=>gs,PickerList:()=>hs,PickerListItem:()=>cs,PickerMenu:()=>as,PickerMenuOption:()=>ls,PropertyStyleSheetBehavior:()=>Xo,Radio:()=>js,RadioGroup:()=>Ns,Registration:()=>Y,ResolverBuilder:()=>u,ResolverImpl:()=>M,Search:()=>Qs,Select:()=>eo,SelectPosition:()=>We,ServiceLocator:()=>x,Skeleton:()=>no,Slider:()=>vo,SliderLabel:()=>ho,SliderMode:()=>mo,StartEnd:()=>o,Switch:()=>yo,Tab:()=>$o,TabPanel:()=>xo,Tabs:()=>Eo,TabsOrientation:()=>ko,TextArea:()=>So,TextAreaResize:()=>To,TextField:()=>Ss,TextFieldType:()=>Rs,Toolbar:()=>Po,Tooltip:()=>Vo,TooltipPosition:()=>zo,TreeItem:()=>qo,TreeView:()=>jo,accordionItemTemplate:()=>h,accordionTemplate:()=>ut,all:()=>O,anchorTemplate:()=>Dt,anchoredRegionTemplate:()=>Ft,applyMixins:()=>dt,avatarTemplate:()=>Xt,badgeTemplate:()=>Zt,breadcrumbItemTemplate:()=>te,breadcrumbTemplate:()=>ie,buttonTemplate:()=>oe,calendarCellTemplate:()=>Te,calendarRowTemplate:()=>Oe,calendarTemplate:()=>Se,calendarWeekdayTemplate:()=>Ee,cardTemplate:()=>Ae,checkboxTemplate:()=>Le,comboboxTemplate:()=>Je,composedContains:()=>oi,composedParent:()=>si,darkModeStylesheetBehavior:()=>Go,dataGridCellTemplate:()=>ii,dataGridRowTemplate:()=>ei,dataGridTemplate:()=>ti,dialogTemplate:()=>Ri,disabledCursor:()=>Qo,disclosureTemplate:()=>Gi,display:()=>Jo,dividerTemplate:()=>Xi,endSlotTemplate:()=>n,endTemplate:()=>r,flipperTemplate:()=>es,focusVisible:()=>tn,forcedColorsStylesheetBehavior:()=>Wo,getDirection:()=>Bt,hidden:()=>Zo,horizontalScrollTemplate:()=>Ks,ignore:()=>A,inject:()=>$,interactiveCalendarGridTemplate:()=>De,isListboxOption:()=>qe,isTreeItemElement:()=>Bo,lazy:()=>R,lightModeStylesheetBehavior:()=>Yo,listboxOptionTemplate:()=>ss,listboxTemplate:()=>ns,menuItemTemplate:()=>$s,menuTemplate:()=>ks,newInstanceForScope:()=>F,newInstanceOf:()=>L,noninteractiveCalendarTemplate:()=>Re,numberFieldTemplate:()=>Ts,optional:()=>S,pickerListItemTemplate:()=>Cs,pickerListTemplate:()=>ys,pickerMenuOptionTemplate:()=>fs,pickerMenuTemplate:()=>bs,pickerTemplate:()=>us,progressRingTemplate:()=>Ms,progressTemplate:()=>zs,radioGroupTemplate:()=>Vs,radioTemplate:()=>Bs,reflectAttributes:()=>Wi,roleForMenuItem:()=>ws,searchTemplate:()=>Gs,selectTemplate:()=>so,singleton:()=>T,skeletonTemplate:()=>oo,sliderLabelTemplate:()=>ao,sliderTemplate:()=>co,startSlotTemplate:()=>a,startTemplate:()=>l,supportsElementInternals:()=>re,switchTemplate:()=>go,tabPanelTemplate:()=>Co,tabTemplate:()=>wo,tabsTemplate:()=>Io,textAreaTemplate:()=>Oo,textFieldTemplate:()=>Ao,toolbarTemplate:()=>Fo,tooltipTemplate:()=>Ho,transient:()=>k,treeItemTemplate:()=>No,treeViewTemplate:()=>Uo,validateKey:()=>X,whitespaceFilter:()=>Ws});var s=i(1351);class o{handleStartContentChange(){this.startContainer.classList.toggle("start",this.start.assignedNodes().length>0)}handleEndContentChange(){this.endContainer.classList.toggle("end",this.end.assignedNodes().length>0)}}const n=(t,e)=>s.html`
    <span
        part="end"
        ${(0,s.ref)("endContainer")}
        class=${t=>e.end?"end":void 0}
    >
        <slot name="end" ${(0,s.ref)("end")} @slotchange="${t=>t.handleEndContentChange()}">
            ${e.end||""}
        </slot>
    </span>
`,a=(t,e)=>s.html`
    <span
        part="start"
        ${(0,s.ref)("startContainer")}
        class="${t=>e.start?"start":void 0}"
    >
        <slot
            name="start"
            ${(0,s.ref)("start")}
            @slotchange="${t=>t.handleStartContentChange()}"
        >
            ${e.start||""}
        </slot>
    </span>
`,r=s.html`
    <span part="end" ${(0,s.ref)("endContainer")}>
        <slot
            name="end"
            ${(0,s.ref)("end")}
            @slotchange="${t=>t.handleEndContentChange()}"
        ></slot>
    </span>
`,l=s.html`
    <span part="start" ${(0,s.ref)("startContainer")}>
        <slot
            name="start"
            ${(0,s.ref)("start")}
            @slotchange="${t=>t.handleStartContentChange()}"
        ></slot>
    </span>
`,h=(t,e)=>s.html`
    <template class="${t=>t.expanded?"expanded":""}">
        <div
            class="heading"
            part="heading"
            role="heading"
            aria-level="${t=>t.headinglevel}"
        >
            <button
                class="button"
                part="button"
                ${(0,s.ref)("expandbutton")}
                aria-expanded="${t=>t.expanded}"
                aria-controls="${t=>t.id}-panel"
                id="${t=>t.id}"
                @click="${(t,e)=>t.clickHandler(e.event)}"
            >
                <span class="heading-content" part="heading-content">
                    <slot name="heading"></slot>
                </span>
            </button>
            ${a(t,e)}
            ${n(t,e)}
            <span class="icon" part="icon" aria-hidden="true">
                <slot name="expanded-icon" part="expanded-icon">
                    ${e.expandedIcon||""}
                </slot>
                <slot name="collapsed-icon" part="collapsed-icon">
                    ${e.collapsedIcon||""}
                </slot>
            <span>
        </div>
        <div
            class="region"
            part="region"
            id="${t=>t.id}-panel"
            role="region"
            aria-labelledby="${t=>t.id}"
        >
            <slot></slot>
        </div>
    </template>
`;var d=i(655);const c=new Map;"metadata"in Reflect||(Reflect.metadata=function(t,e){return function(i){Reflect.defineMetadata(t,e,i)}},Reflect.defineMetadata=function(t,e,i){let s=c.get(i);void 0===s&&c.set(i,s=new Map),s.set(t,e)},Reflect.getOwnMetadata=function(t,e){const i=c.get(e);if(void 0!==i)return i.get(t)});class u{constructor(t,e){this.container=t,this.key=e}instance(t){return this.registerResolver(0,t)}singleton(t){return this.registerResolver(1,t)}transient(t){return this.registerResolver(2,t)}callback(t){return this.registerResolver(3,t)}cachedCallback(t){return this.registerResolver(3,G(t))}aliasTo(t){return this.registerResolver(5,t)}registerResolver(t,e){const{container:i,key:s}=this;return this.container=this.key=void 0,i.registerResolver(s,new M(s,t,e))}}function p(t){const e=t.slice(),i=Object.keys(t),s=i.length;let o;for(let n=0;n<s;++n)o=i[n],it(o)||(e[o]=t[o]);return e}const m=Object.freeze({none(t){throw Error(`${t.toString()} not registered, did you forget to add @singleton()?`)},singleton:t=>new M(t,1,t),transient:t=>new M(t,2,t)}),v=Object.freeze({default:Object.freeze({parentLocator:()=>null,responsibleForOwnerRequests:!1,defaultResolver:m.singleton})}),g=new Map;function b(t){return e=>Reflect.getOwnMetadata(t,e)}let f=null;const y=Object.freeze({createContainer:t=>new K(null,Object.assign({},v.default,t)),findResponsibleContainer(t){const e=t.$$container$$;return e&&e.responsibleForOwnerRequests?e:y.findParentContainer(t)},findParentContainer(t){const e=new CustomEvent(j,{bubbles:!0,composed:!0,cancelable:!0,detail:{container:void 0}});return t.dispatchEvent(e),e.detail.container||y.getOrCreateDOMContainer()},getOrCreateDOMContainer:(t,e)=>t?t.$$container$$||new K(t,Object.assign({},v.default,e,{parentLocator:y.findParentContainer})):f||(f=new K(null,Object.assign({},v.default,e,{parentLocator:()=>null}))),getDesignParamtypes:b("design:paramtypes"),getAnnotationParamtypes:b("di:paramtypes"),getOrCreateAnnotationParamTypes(t){let e=this.getAnnotationParamtypes(t);return void 0===e&&Reflect.defineMetadata("di:paramtypes",e=[],t),e},getDependencies(t){let e=g.get(t);if(void 0===e){const i=t.inject;if(void 0===i){const i=y.getDesignParamtypes(t),s=y.getAnnotationParamtypes(t);if(void 0===i)if(void 0===s){const i=Object.getPrototypeOf(t);e="function"==typeof i&&i!==Function.prototype?p(y.getDependencies(i)):[]}else e=p(s);else if(void 0===s)e=p(i);else{e=p(i);let t,o=s.length;for(let i=0;i<o;++i)t=s[i],void 0!==t&&(e[i]=t);const n=Object.keys(s);let a;o=n.length;for(let t=0;t<o;++t)a=n[t],it(a)||(e[a]=s[a])}}else e=p(i);g.set(t,e)}return e},defineProperty(t,e,i,o=!1){const n=`$di_${e}`;Reflect.defineProperty(t,e,{get:function(){let t=this[n];if(void 0===t){const a=this instanceof HTMLElement?y.findResponsibleContainer(this):y.getOrCreateDOMContainer();if(t=a.get(i),this[n]=t,o&&this instanceof s.FASTElement){const s=this.$fastController,o=()=>{y.findResponsibleContainer(this).get(i)!==this[n]&&(this[n]=t,s.notify(e))};s.subscribe({handleChange:o},"isConnected")}}return t}})},createInterface(t,e){const i="function"==typeof t?t:e,s="string"==typeof t?t:t&&"friendlyName"in t&&t.friendlyName||Z,o="string"!=typeof t&&(t&&"respectConnection"in t&&t.respectConnection||!1),n=function(t,e,i){if(null==t||void 0!==new.target)throw new Error(`No registration for interface: '${n.friendlyName}'`);e?y.defineProperty(t,e,n,o):y.getOrCreateAnnotationParamTypes(t)[i]=n};return n.$isInterface=!0,n.friendlyName=null==s?"(anonymous)":s,null!=i&&(n.register=function(t,e){return i(new u(t,null!=e?e:n))}),n.toString=function(){return`InterfaceSymbol<${n.friendlyName}>`},n},inject:(...t)=>function(e,i,s){if("number"==typeof s){const i=y.getOrCreateAnnotationParamTypes(e),o=t[0];void 0!==o&&(i[s]=o)}else if(i)y.defineProperty(e,i,t[0]);else{const i=s?y.getOrCreateAnnotationParamTypes(s.value):y.getOrCreateAnnotationParamTypes(e);let o;for(let e=0;e<t.length;++e)o=t[e],void 0!==o&&(i[e]=o)}},transient:t=>(t.register=function(e){return Y.transient(t,t).register(e)},t.registerInRequestor=!1,t),singleton:(t,e=E)=>(t.register=function(e){return Y.singleton(t,t).register(e)},t.registerInRequestor=e.scoped,t)}),C=y.createInterface("Container"),x=C;function w(t){return function(e){const i=function(t,e,s){y.inject(i)(t,e,s)};return i.$isResolver=!0,i.resolve=function(i,s){return t(e,i,s)},i}}const $=y.inject;function I(t){return y.transient(t)}function k(t){return null==t?I:I(t)}const E={scoped:!1};function T(t){return"function"==typeof t?y.singleton(t):function(e){return y.singleton(e,t)}}const O=(D=(t,e,i,s)=>i.getAll(t,s),function(t,e){e=!!e;const i=function(t,e,s){y.inject(i)(t,e,s)};return i.$isResolver=!0,i.resolve=function(i,s){return D(t,0,s,e)},i});var D;const R=w(((t,e,i)=>()=>i.get(t))),S=w(((t,e,i)=>i.has(t,!0)?i.get(t):void 0));function A(t,e,i){y.inject(A)(t,e,i)}A.$isResolver=!0,A.resolve=()=>{};const F=w(((t,e,i)=>{const s=P(t,e),o=new M(t,0,s);return i.registerResolver(t,o),s})),L=w(((t,e,i)=>P(t,e)));function P(t,e){return e.getFactory(t).construct(e)}class M{constructor(t,e,i){this.key=t,this.strategy=e,this.state=i,this.resolving=!1}get $isResolver(){return!0}register(t){return t.registerResolver(this.key,this)}resolve(t,e){switch(this.strategy){case 0:return this.state;case 1:if(this.resolving)throw new Error(`Cyclic dependency found: ${this.state.name}`);return this.resolving=!0,this.state=t.getFactory(this.state).construct(e),this.strategy=0,this.resolving=!1,this.state;case 2:{const i=t.getFactory(this.state);if(null===i)throw new Error(`Resolver for ${String(this.key)} returned a null factory`);return i.construct(e)}case 3:return this.state(t,e,this);case 4:return this.state[0].resolve(t,e);case 5:return e.get(this.state);default:throw new Error(`Invalid resolver strategy specified: ${this.strategy}.`)}}getFactory(t){var e,i,s;switch(this.strategy){case 1:case 2:return t.getFactory(this.state);case 5:return null!==(s=null===(i=null===(e=t.getResolver(this.state))||void 0===e?void 0:e.getFactory)||void 0===i?void 0:i.call(e,t))&&void 0!==s?s:null;default:return null}}}function H(t){return this.get(t)}function z(t,e){return e(t)}class V{constructor(t,e){this.Type=t,this.dependencies=e,this.transformers=null}construct(t,e){let i;return i=void 0===e?new this.Type(...this.dependencies.map(H,t)):new this.Type(...this.dependencies.map(H,t),...e),null==this.transformers?i:this.transformers.reduce(z,i)}registerTransformer(t){(this.transformers||(this.transformers=[])).push(t)}}const N={$isResolver:!0,resolve:(t,e)=>e};function B(t){return"function"==typeof t.register}function q(t){return function(t){return B(t)&&"boolean"==typeof t.registerInRequestor}(t)&&t.registerInRequestor}const U=new Set(["Array","ArrayBuffer","Boolean","DataView","Date","Error","EvalError","Float32Array","Float64Array","Function","Int8Array","Int16Array","Int32Array","Map","Number","Object","Promise","RangeError","ReferenceError","RegExp","Set","SharedArrayBuffer","String","SyntaxError","TypeError","Uint8Array","Uint8ClampedArray","Uint16Array","Uint32Array","URIError","WeakMap","WeakSet"]),j="__DI_LOCATE_PARENT__",_=new Map;class K{constructor(t,e){this.owner=t,this.config=e,this._parent=void 0,this.registerDepth=0,this.context=null,null!==t&&(t.$$container$$=this),this.resolvers=new Map,this.resolvers.set(C,N),t instanceof Node&&t.addEventListener(j,(t=>{t.composedPath()[0]!==this.owner&&(t.detail.container=this,t.stopImmediatePropagation())}))}get parent(){return void 0===this._parent&&(this._parent=this.config.parentLocator(this.owner)),this._parent}get depth(){return null===this.parent?0:this.parent.depth+1}get responsibleForOwnerRequests(){return this.config.responsibleForOwnerRequests}registerWithContext(t,...e){return this.context=t,this.register(...e),this.context=null,this}register(...t){if(100==++this.registerDepth)throw new Error("Unable to autoregister dependency");let e,i,s,o,n;const a=this.context;for(let r=0,l=t.length;r<l;++r)if(e=t[r],J(e))if(B(e))e.register(this,a);else if(void 0!==e.prototype)Y.singleton(e,e).register(this);else for(i=Object.keys(e),o=0,n=i.length;o<n;++o)s=e[i[o]],J(s)&&(B(s)?s.register(this,a):this.register(s));return--this.registerDepth,this}registerResolver(t,e){X(t);const i=this.resolvers,s=i.get(t);return null==s?i.set(t,e):s instanceof M&&4===s.strategy?s.state.push(e):i.set(t,new M(t,4,[s,e])),e}registerTransformer(t,e){const i=this.getResolver(t);if(null==i)return!1;if(i.getFactory){const t=i.getFactory(this);return null!=t&&(t.registerTransformer(e),!0)}return!1}getResolver(t,e=!0){if(X(t),void 0!==t.resolve)return t;let i,s=this;for(;null!=s;){if(i=s.resolvers.get(t),null!=i)return i;if(null==s.parent){const i=q(t)?this:s;return e?this.jitRegister(t,i):null}s=s.parent}return null}has(t,e=!1){return!!this.resolvers.has(t)||!(!e||null==this.parent)&&this.parent.has(t,!0)}get(t){if(X(t),t.$isResolver)return t.resolve(this,this);let e,i=this;for(;null!=i;){if(e=i.resolvers.get(t),null!=e)return e.resolve(i,this);if(null==i.parent){const s=q(t)?this:i;return e=this.jitRegister(t,s),e.resolve(i,this)}i=i.parent}throw new Error(`Unable to resolve key: ${String(t)}`)}getAll(t,e=!1){X(t);const i=this;let o,n=i;if(e){let e=s.emptyArray;for(;null!=n;)o=n.resolvers.get(t),null!=o&&(e=e.concat(Q(o,n,i))),n=n.parent;return e}for(;null!=n;){if(o=n.resolvers.get(t),null!=o)return Q(o,n,i);if(n=n.parent,null==n)return s.emptyArray}return s.emptyArray}getFactory(t){let e=_.get(t);if(void 0===e){if(tt(t))throw new Error(`${t.name} is a native function and therefore cannot be safely constructed by DI. If this is intentional, please use a callback or cachedCallback resolver.`);_.set(t,e=new V(t,y.getDependencies(t)))}return e}registerFactory(t,e){_.set(t,e)}createChild(t){return new K(null,Object.assign({},this.config,t,{parentLocator:()=>this}))}jitRegister(t,e){if("function"!=typeof t)throw new Error(`Attempted to jitRegister something that is not a constructor: '${t}'. Did you forget to register this dependency?`);if(U.has(t.name))throw new Error(`Attempted to jitRegister an intrinsic type: ${t.name}. Did you forget to add @inject(Key)`);if(B(t)){const i=t.register(e);if(!(i instanceof Object)||null==i.resolve){const i=e.resolvers.get(t);if(null!=i)return i;throw new Error("A valid resolver was not returned from the static register method")}return i}if(t.$isInterface)throw new Error(`Attempted to jitRegister an interface: ${t.friendlyName}`);{const i=this.config.defaultResolver(t,e);return e.resolvers.set(t,i),i}}}const W=new WeakMap;function G(t){return function(e,i,s){if(W.has(s))return W.get(s);const o=t(e,i,s);return W.set(s,o),o}}const Y=Object.freeze({instance:(t,e)=>new M(t,0,e),singleton:(t,e)=>new M(t,1,e),transient:(t,e)=>new M(t,2,e),callback:(t,e)=>new M(t,3,e),cachedCallback:(t,e)=>new M(t,3,G(e)),aliasTo:(t,e)=>new M(e,5,t)});function X(t){if(null==t)throw new Error("key/value cannot be null or undefined. Are you trying to inject/register something that doesn't exist with DI?")}function Q(t,e,i){if(t instanceof M&&4===t.strategy){const s=t.state;let o=s.length;const n=new Array(o);for(;o--;)n[o]=s[o].resolve(e,i);return n}return[t.resolve(e,i)]}const Z="(anonymous)";function J(t){return"object"==typeof t&&null!==t||"function"==typeof t}const tt=function(){const t=new WeakMap;let e=!1,i="",s=0;return function(o){return e=t.get(o),void 0===e&&(i=o.toString(),s=i.length,e=s>=29&&s<=100&&125===i.charCodeAt(s-1)&&i.charCodeAt(s-2)<=32&&93===i.charCodeAt(s-3)&&101===i.charCodeAt(s-4)&&100===i.charCodeAt(s-5)&&111===i.charCodeAt(s-6)&&99===i.charCodeAt(s-7)&&32===i.charCodeAt(s-8)&&101===i.charCodeAt(s-9)&&118===i.charCodeAt(s-10)&&105===i.charCodeAt(s-11)&&116===i.charCodeAt(s-12)&&97===i.charCodeAt(s-13)&&110===i.charCodeAt(s-14)&&88===i.charCodeAt(s-15),t.set(o,e)),e}}(),et={};function it(t){switch(typeof t){case"number":return t>=0&&(0|t)===t;case"string":{const e=et[t];if(void 0!==e)return e;const i=t.length;if(0===i)return et[t]=!1;let s=0;for(let e=0;e<i;++e)if(s=t.charCodeAt(e),0===e&&48===s&&i>1||s<48||s>57)return et[t]=!1;return et[t]=!0}default:return!1}}function st(t){return`${t.toLowerCase()}:presentation`}const ot=new Map,nt=Object.freeze({define(t,e,i){const s=st(t);void 0===ot.get(s)?ot.set(s,e):ot.set(s,!1),i.register(Y.instance(s,e))},forTag(t,e){const i=st(t),s=ot.get(i);return!1===s?y.findResponsibleContainer(e).get(i):s||null}});class at{constructor(t,e){this.template=t||null,this.styles=void 0===e?null:Array.isArray(e)?s.ElementStyles.create(e):e instanceof s.ElementStyles?e:s.ElementStyles.create([e])}applyTo(t){const e=t.$fastController;null===e.template&&(e.template=this.template),null===e.styles&&(e.styles=this.styles)}}class rt extends s.FASTElement{constructor(){super(...arguments),this._presentation=void 0}get $presentation(){return void 0===this._presentation&&(this._presentation=nt.forTag(this.tagName,this)),this._presentation}templateChanged(){void 0!==this.template&&(this.$fastController.template=this.template)}stylesChanged(){void 0!==this.styles&&(this.$fastController.styles=this.styles)}connectedCallback(){null!==this.$presentation&&this.$presentation.applyTo(this),super.connectedCallback()}static compose(t){return(e={})=>new ht(this===rt?class extends rt{}:this,t,e)}}function lt(t,e,i){return"function"==typeof t?t(e,i):t}(0,d.gn)([s.observable],rt.prototype,"template",void 0),(0,d.gn)([s.observable],rt.prototype,"styles",void 0);class ht{constructor(t,e,i){this.type=t,this.elementDefinition=e,this.overrideDefinition=i,this.definition=Object.assign(Object.assign({},this.elementDefinition),this.overrideDefinition)}register(t,e){const i=this.definition,s=this.overrideDefinition,o=`${i.prefix||e.elementPrefix}-${i.baseName}`;e.tryDefineElement({name:o,type:this.type,baseClass:this.elementDefinition.baseClass,callback:t=>{const e=new at(lt(i.template,t,i),lt(i.styles,t,i));t.definePresentation(e);let o=lt(i.shadowOptions,t,i);t.shadowRootMode&&(o?s.shadowOptions||(o.mode=t.shadowRootMode):null!==o&&(o={mode:t.shadowRootMode})),t.defineElement({elementOptions:lt(i.elementOptions,t,i),shadowOptions:o,attributes:lt(i.attributes,t,i)})}})}}function dt(t,...e){const i=s.AttributeConfiguration.locate(t);e.forEach((e=>{Object.getOwnPropertyNames(e.prototype).forEach((i=>{"constructor"!==i&&Object.defineProperty(t.prototype,i,Object.getOwnPropertyDescriptor(e.prototype,i))})),s.AttributeConfiguration.locate(e).forEach((t=>i.push(t)))}))}class ct extends rt{constructor(){super(...arguments),this.headinglevel=2,this.expanded=!1,this.clickHandler=t=>{this.expanded=!this.expanded,this.change()},this.change=()=>{this.$emit("change")}}}(0,d.gn)([(0,s.attr)({attribute:"heading-level",mode:"fromView",converter:s.nullableNumberConverter})],ct.prototype,"headinglevel",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],ct.prototype,"expanded",void 0),(0,d.gn)([s.attr],ct.prototype,"id",void 0),dt(ct,o);const ut=(t,e)=>s.html`
    <template>
        <slot ${(0,s.slotted)({property:"accordionItems",filter:(0,s.elements)()})}></slot>
        <slot name="item" part="item" ${(0,s.slotted)("accordionItems")}></slot>
    </template>
`;var pt;!function(t){t[t.alt=18]="alt",t[t.arrowDown=40]="arrowDown",t[t.arrowLeft=37]="arrowLeft",t[t.arrowRight=39]="arrowRight",t[t.arrowUp=38]="arrowUp",t[t.back=8]="back",t[t.backSlash=220]="backSlash",t[t.break=19]="break",t[t.capsLock=20]="capsLock",t[t.closeBracket=221]="closeBracket",t[t.colon=186]="colon",t[t.colon2=59]="colon2",t[t.comma=188]="comma",t[t.ctrl=17]="ctrl",t[t.delete=46]="delete",t[t.end=35]="end",t[t.enter=13]="enter",t[t.equals=187]="equals",t[t.equals2=61]="equals2",t[t.equals3=107]="equals3",t[t.escape=27]="escape",t[t.forwardSlash=191]="forwardSlash",t[t.function1=112]="function1",t[t.function10=121]="function10",t[t.function11=122]="function11",t[t.function12=123]="function12",t[t.function2=113]="function2",t[t.function3=114]="function3",t[t.function4=115]="function4",t[t.function5=116]="function5",t[t.function6=117]="function6",t[t.function7=118]="function7",t[t.function8=119]="function8",t[t.function9=120]="function9",t[t.home=36]="home",t[t.insert=45]="insert",t[t.menu=93]="menu",t[t.minus=189]="minus",t[t.minus2=109]="minus2",t[t.numLock=144]="numLock",t[t.numPad0=96]="numPad0",t[t.numPad1=97]="numPad1",t[t.numPad2=98]="numPad2",t[t.numPad3=99]="numPad3",t[t.numPad4=100]="numPad4",t[t.numPad5=101]="numPad5",t[t.numPad6=102]="numPad6",t[t.numPad7=103]="numPad7",t[t.numPad8=104]="numPad8",t[t.numPad9=105]="numPad9",t[t.numPadDivide=111]="numPadDivide",t[t.numPadDot=110]="numPadDot",t[t.numPadMinus=109]="numPadMinus",t[t.numPadMultiply=106]="numPadMultiply",t[t.numPadPlus=107]="numPadPlus",t[t.openBracket=219]="openBracket",t[t.pageDown=34]="pageDown",t[t.pageUp=33]="pageUp",t[t.period=190]="period",t[t.print=44]="print",t[t.quote=222]="quote",t[t.scrollLock=145]="scrollLock",t[t.shift=16]="shift",t[t.space=32]="space",t[t.tab=9]="tab",t[t.tilde=192]="tilde",t[t.windowsLeft=91]="windowsLeft",t[t.windowsOpera=219]="windowsOpera",t[t.windowsRight=92]="windowsRight"}(pt||(pt={}));const mt="ArrowDown",vt="ArrowLeft",gt="ArrowRight",bt="ArrowUp",ft="Enter",yt="Escape",Ct="Home",xt="End",wt=" ",$t="Tab",It={ArrowDown:mt,ArrowLeft:vt,ArrowRight:gt,ArrowUp:bt};function kt(t,e,i){return Math.min(Math.max(i,t),e)}function Et(t,e,i=0){return[e,i]=[e,i].sort(((t,e)=>t-e)),e<=t&&t<i}const Tt={single:"single",multi:"multi"};class Ot extends rt{constructor(){super(...arguments),this.expandmode=Tt.multi,this.activeItemIndex=0,this.change=()=>{this.$emit("change",this.activeid)},this.setItems=()=>{var t;0!==this.accordionItems.length&&(this.accordionIds=this.getItemIds(),this.accordionItems.forEach(((t,e)=>{t instanceof ct&&(t.addEventListener("change",this.activeItemChange),this.isSingleExpandMode()&&(this.activeItemIndex!==e?t.expanded=!1:t.expanded=!0));const i=this.accordionIds[e];t.setAttribute("id","string"!=typeof i?`accordion-${e+1}`:i),this.activeid=this.accordionIds[this.activeItemIndex],t.addEventListener("keydown",this.handleItemKeyDown),t.addEventListener("focus",this.handleItemFocus)})),this.isSingleExpandMode())&&(null!==(t=this.findExpandedItem())&&void 0!==t?t:this.accordionItems[0]).setAttribute("aria-disabled","true")},this.removeItemListeners=t=>{t.forEach(((t,e)=>{t.removeEventListener("change",this.activeItemChange),t.removeEventListener("keydown",this.handleItemKeyDown),t.removeEventListener("focus",this.handleItemFocus)}))},this.activeItemChange=t=>{if(t.defaultPrevented||t.target!==t.currentTarget)return;t.preventDefault();const e=t.target;this.activeid=e.getAttribute("id"),this.isSingleExpandMode()&&(this.resetItems(),e.expanded=!0,e.setAttribute("aria-disabled","true"),this.accordionItems.forEach((t=>{t.hasAttribute("disabled")||t.id===this.activeid||t.removeAttribute("aria-disabled")}))),this.activeItemIndex=Array.from(this.accordionItems).indexOf(e),this.change()},this.handleItemKeyDown=t=>{if(t.target===t.currentTarget)switch(this.accordionIds=this.getItemIds(),t.key){case bt:t.preventDefault(),this.adjust(-1);break;case mt:t.preventDefault(),this.adjust(1);break;case Ct:this.activeItemIndex=0,this.focusItem();break;case xt:this.activeItemIndex=this.accordionItems.length-1,this.focusItem()}},this.handleItemFocus=t=>{if(t.target===t.currentTarget){const e=t.target,i=this.activeItemIndex=Array.from(this.accordionItems).indexOf(e);this.activeItemIndex!==i&&-1!==i&&(this.activeItemIndex=i,this.activeid=this.accordionIds[this.activeItemIndex])}}}accordionItemsChanged(t,e){this.$fastController.isConnected&&(this.removeItemListeners(t),this.setItems())}findExpandedItem(){for(let t=0;t<this.accordionItems.length;t++)if("true"===this.accordionItems[t].getAttribute("expanded"))return this.accordionItems[t];return null}resetItems(){this.accordionItems.forEach(((t,e)=>{t.expanded=!1}))}getItemIds(){return this.accordionItems.map((t=>t.getAttribute("id")))}isSingleExpandMode(){return this.expandmode===Tt.single}adjust(t){var e,i;this.activeItemIndex=(0,e=this.accordionItems.length-1,(i=this.activeItemIndex+t)<0?e:i>e?0:i),this.focusItem()}focusItem(){const t=this.accordionItems[this.activeItemIndex];t instanceof ct&&t.expandbutton.focus()}}(0,d.gn)([(0,s.attr)({attribute:"expand-mode"})],Ot.prototype,"expandmode",void 0),(0,d.gn)([s.observable],Ot.prototype,"accordionItems",void 0);const Dt=(t,e)=>s.html`
    <a
        class="control"
        part="control"
        download="${t=>t.download}"
        href="${t=>t.href}"
        hreflang="${t=>t.hreflang}"
        ping="${t=>t.ping}"
        referrerpolicy="${t=>t.referrerpolicy}"
        rel="${t=>t.rel}"
        target="${t=>t.target}"
        type="${t=>t.type}"
        aria-atomic="${t=>t.ariaAtomic}"
        aria-busy="${t=>t.ariaBusy}"
        aria-controls="${t=>t.ariaControls}"
        aria-current="${t=>t.ariaCurrent}"
        aria-describedby="${t=>t.ariaDescribedby}"
        aria-details="${t=>t.ariaDetails}"
        aria-disabled="${t=>t.ariaDisabled}"
        aria-errormessage="${t=>t.ariaErrormessage}"
        aria-expanded="${t=>t.ariaExpanded}"
        aria-flowto="${t=>t.ariaFlowto}"
        aria-haspopup="${t=>t.ariaHaspopup}"
        aria-hidden="${t=>t.ariaHidden}"
        aria-invalid="${t=>t.ariaInvalid}"
        aria-keyshortcuts="${t=>t.ariaKeyshortcuts}"
        aria-label="${t=>t.ariaLabel}"
        aria-labelledby="${t=>t.ariaLabelledby}"
        aria-live="${t=>t.ariaLive}"
        aria-owns="${t=>t.ariaOwns}"
        aria-relevant="${t=>t.ariaRelevant}"
        aria-roledescription="${t=>t.ariaRoledescription}"
        ${(0,s.ref)("control")}
    >
        ${a(t,e)}
        <span class="content" part="content">
            <slot ${(0,s.slotted)("defaultSlottedContent")}></slot>
        </span>
        ${n(t,e)}
    </a>
`;class Rt{}(0,d.gn)([(0,s.attr)({attribute:"aria-atomic"})],Rt.prototype,"ariaAtomic",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-busy"})],Rt.prototype,"ariaBusy",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-controls"})],Rt.prototype,"ariaControls",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-current"})],Rt.prototype,"ariaCurrent",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-describedby"})],Rt.prototype,"ariaDescribedby",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-details"})],Rt.prototype,"ariaDetails",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-disabled"})],Rt.prototype,"ariaDisabled",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-errormessage"})],Rt.prototype,"ariaErrormessage",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-flowto"})],Rt.prototype,"ariaFlowto",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-haspopup"})],Rt.prototype,"ariaHaspopup",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-hidden"})],Rt.prototype,"ariaHidden",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-invalid"})],Rt.prototype,"ariaInvalid",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-keyshortcuts"})],Rt.prototype,"ariaKeyshortcuts",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-label"})],Rt.prototype,"ariaLabel",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-labelledby"})],Rt.prototype,"ariaLabelledby",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-live"})],Rt.prototype,"ariaLive",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-owns"})],Rt.prototype,"ariaOwns",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-relevant"})],Rt.prototype,"ariaRelevant",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-roledescription"})],Rt.prototype,"ariaRoledescription",void 0);class St extends rt{constructor(){super(...arguments),this.handleUnsupportedDelegatesFocus=()=>{var t;window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&(null===(t=this.$fastController.definition.shadowOptions)||void 0===t?void 0:t.delegatesFocus)&&(this.focus=()=>{var t;null===(t=this.control)||void 0===t||t.focus()})}}connectedCallback(){super.connectedCallback(),this.handleUnsupportedDelegatesFocus()}}(0,d.gn)([s.attr],St.prototype,"download",void 0),(0,d.gn)([s.attr],St.prototype,"href",void 0),(0,d.gn)([s.attr],St.prototype,"hreflang",void 0),(0,d.gn)([s.attr],St.prototype,"ping",void 0),(0,d.gn)([s.attr],St.prototype,"referrerpolicy",void 0),(0,d.gn)([s.attr],St.prototype,"rel",void 0),(0,d.gn)([s.attr],St.prototype,"target",void 0),(0,d.gn)([s.attr],St.prototype,"type",void 0),(0,d.gn)([s.observable],St.prototype,"defaultSlottedContent",void 0);class At{}(0,d.gn)([(0,s.attr)({attribute:"aria-expanded"})],At.prototype,"ariaExpanded",void 0),dt(At,Rt),dt(St,o,At);const Ft=(t,e)=>s.html`
    <template class="${t=>t.initialLayoutComplete?"loaded":""}">
        ${(0,s.when)((t=>t.initialLayoutComplete),s.html`
                <slot></slot>
            `)}
    </template>
`;var Lt=i(2586);const Pt="focus",Mt="focusin",Ht="focusout",zt="keydown",Vt="resize",Nt="scroll",Bt=t=>{const e=t.closest("[dir]");return null!==e&&"rtl"===e.dir?Lt.N.rtl:Lt.N.ltr};class qt extends rt{constructor(){super(...arguments),this.anchor="",this.viewport="",this.horizontalPositioningMode="uncontrolled",this.horizontalDefaultPosition="unset",this.horizontalViewportLock=!1,this.horizontalInset=!1,this.horizontalScaling="content",this.verticalPositioningMode="uncontrolled",this.verticalDefaultPosition="unset",this.verticalViewportLock=!1,this.verticalInset=!1,this.verticalScaling="content",this.fixedPlacement=!1,this.autoUpdateMode="anchor",this.anchorElement=null,this.viewportElement=null,this.initialLayoutComplete=!1,this.resizeDetector=null,this.baseHorizontalOffset=0,this.baseVerticalOffset=0,this.pendingPositioningUpdate=!1,this.pendingReset=!1,this.currentDirection=Lt.N.ltr,this.regionVisible=!1,this.forceUpdate=!1,this.updateThreshold=.5,this.update=()=>{this.pendingPositioningUpdate||this.requestPositionUpdates()},this.startObservers=()=>{this.stopObservers(),null!==this.anchorElement&&(this.requestPositionUpdates(),null!==this.resizeDetector&&(this.resizeDetector.observe(this.anchorElement),this.resizeDetector.observe(this)))},this.requestPositionUpdates=()=>{null===this.anchorElement||this.pendingPositioningUpdate||(qt.intersectionService.requestPosition(this,this.handleIntersection),qt.intersectionService.requestPosition(this.anchorElement,this.handleIntersection),null!==this.viewportElement&&qt.intersectionService.requestPosition(this.viewportElement,this.handleIntersection),this.pendingPositioningUpdate=!0)},this.stopObservers=()=>{this.pendingPositioningUpdate&&(this.pendingPositioningUpdate=!1,qt.intersectionService.cancelRequestPosition(this,this.handleIntersection),null!==this.anchorElement&&qt.intersectionService.cancelRequestPosition(this.anchorElement,this.handleIntersection),null!==this.viewportElement&&qt.intersectionService.cancelRequestPosition(this.viewportElement,this.handleIntersection)),null!==this.resizeDetector&&this.resizeDetector.disconnect()},this.getViewport=()=>"string"!=typeof this.viewport||""===this.viewport?document.documentElement:document.getElementById(this.viewport),this.getAnchor=()=>document.getElementById(this.anchor),this.handleIntersection=t=>{this.pendingPositioningUpdate&&(this.pendingPositioningUpdate=!1,this.applyIntersectionEntries(t)&&this.updateLayout())},this.applyIntersectionEntries=t=>{const e=t.find((t=>t.target===this)),i=t.find((t=>t.target===this.anchorElement)),s=t.find((t=>t.target===this.viewportElement));return void 0!==e&&void 0!==s&&void 0!==i&&!!(!this.regionVisible||this.forceUpdate||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect||this.isRectDifferent(this.anchorRect,i.boundingClientRect)||this.isRectDifferent(this.viewportRect,s.boundingClientRect)||this.isRectDifferent(this.regionRect,e.boundingClientRect))&&(this.regionRect=e.boundingClientRect,this.anchorRect=i.boundingClientRect,this.viewportElement===document.documentElement?this.viewportRect=new DOMRectReadOnly(s.boundingClientRect.x+document.documentElement.scrollLeft,s.boundingClientRect.y+document.documentElement.scrollTop,s.boundingClientRect.width,s.boundingClientRect.height):this.viewportRect=s.boundingClientRect,this.updateRegionOffset(),this.forceUpdate=!1,!0)},this.updateRegionOffset=()=>{this.anchorRect&&this.regionRect&&(this.baseHorizontalOffset=this.baseHorizontalOffset+(this.anchorRect.left-this.regionRect.left)+(this.translateX-this.baseHorizontalOffset),this.baseVerticalOffset=this.baseVerticalOffset+(this.anchorRect.top-this.regionRect.top)+(this.translateY-this.baseVerticalOffset))},this.isRectDifferent=(t,e)=>Math.abs(t.top-e.top)>this.updateThreshold||Math.abs(t.right-e.right)>this.updateThreshold||Math.abs(t.bottom-e.bottom)>this.updateThreshold||Math.abs(t.left-e.left)>this.updateThreshold,this.handleResize=t=>{this.update()},this.reset=()=>{this.pendingReset&&(this.pendingReset=!1,null===this.anchorElement&&(this.anchorElement=this.getAnchor()),null===this.viewportElement&&(this.viewportElement=this.getViewport()),this.currentDirection=Bt(this),this.startObservers())},this.updateLayout=()=>{let t,e;if("uncontrolled"!==this.horizontalPositioningMode){const t=this.getPositioningOptions(this.horizontalInset);if("center"===this.horizontalDefaultPosition)e="center";else if("unset"!==this.horizontalDefaultPosition){let t=this.horizontalDefaultPosition;if("start"===t||"end"===t){const e=Bt(this);if(e!==this.currentDirection)return this.currentDirection=e,void this.initialize();t=this.currentDirection===Lt.N.ltr?"start"===t?"left":"right":"start"===t?"right":"left"}switch(t){case"left":e=this.horizontalInset?"insetStart":"start";break;case"right":e=this.horizontalInset?"insetEnd":"end"}}const i=void 0!==this.horizontalThreshold?this.horizontalThreshold:void 0!==this.regionRect?this.regionRect.width:0,s=void 0!==this.anchorRect?this.anchorRect.left:0,o=void 0!==this.anchorRect?this.anchorRect.right:0,n=void 0!==this.anchorRect?this.anchorRect.width:0,a=void 0!==this.viewportRect?this.viewportRect.left:0,r=void 0!==this.viewportRect?this.viewportRect.right:0;(void 0===e||"locktodefault"!==this.horizontalPositioningMode&&this.getAvailableSpace(e,s,o,n,a,r)<i)&&(e=this.getAvailableSpace(t[0],s,o,n,a,r)>this.getAvailableSpace(t[1],s,o,n,a,r)?t[0]:t[1])}if("uncontrolled"!==this.verticalPositioningMode){const e=this.getPositioningOptions(this.verticalInset);if("center"===this.verticalDefaultPosition)t="center";else if("unset"!==this.verticalDefaultPosition)switch(this.verticalDefaultPosition){case"top":t=this.verticalInset?"insetStart":"start";break;case"bottom":t=this.verticalInset?"insetEnd":"end"}const i=void 0!==this.verticalThreshold?this.verticalThreshold:void 0!==this.regionRect?this.regionRect.height:0,s=void 0!==this.anchorRect?this.anchorRect.top:0,o=void 0!==this.anchorRect?this.anchorRect.bottom:0,n=void 0!==this.anchorRect?this.anchorRect.height:0,a=void 0!==this.viewportRect?this.viewportRect.top:0,r=void 0!==this.viewportRect?this.viewportRect.bottom:0;(void 0===t||"locktodefault"!==this.verticalPositioningMode&&this.getAvailableSpace(t,s,o,n,a,r)<i)&&(t=this.getAvailableSpace(e[0],s,o,n,a,r)>this.getAvailableSpace(e[1],s,o,n,a,r)?e[0]:e[1])}const i=this.getNextRegionDimension(e,t),s=this.horizontalPosition!==e||this.verticalPosition!==t;if(this.setHorizontalPosition(e,i),this.setVerticalPosition(t,i),this.updateRegionStyle(),!this.initialLayoutComplete)return this.initialLayoutComplete=!0,void this.requestPositionUpdates();this.regionVisible||(this.regionVisible=!0,this.style.removeProperty("pointer-events"),this.style.removeProperty("opacity"),this.classList.toggle("loaded",!0),this.$emit("loaded",this,{bubbles:!1})),this.updatePositionClasses(),s&&this.$emit("positionchange",this,{bubbles:!1})},this.updateRegionStyle=()=>{this.style.width=this.regionWidth,this.style.height=this.regionHeight,this.style.transform=`translate(${this.translateX}px, ${this.translateY}px)`},this.updatePositionClasses=()=>{this.classList.toggle("top","start"===this.verticalPosition),this.classList.toggle("bottom","end"===this.verticalPosition),this.classList.toggle("inset-top","insetStart"===this.verticalPosition),this.classList.toggle("inset-bottom","insetEnd"===this.verticalPosition),this.classList.toggle("vertical-center","center"===this.verticalPosition),this.classList.toggle("left","start"===this.horizontalPosition),this.classList.toggle("right","end"===this.horizontalPosition),this.classList.toggle("inset-left","insetStart"===this.horizontalPosition),this.classList.toggle("inset-right","insetEnd"===this.horizontalPosition),this.classList.toggle("horizontal-center","center"===this.horizontalPosition)},this.setHorizontalPosition=(t,e)=>{if(void 0===t||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect)return;let i=0;switch(this.horizontalScaling){case"anchor":case"fill":i=this.horizontalViewportLock?this.viewportRect.width:e.width,this.regionWidth=`${i}px`;break;case"content":i=this.regionRect.width,this.regionWidth="unset"}let s=0;switch(t){case"start":this.translateX=this.baseHorizontalOffset-i,this.horizontalViewportLock&&this.anchorRect.left>this.viewportRect.right&&(this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.right));break;case"insetStart":this.translateX=this.baseHorizontalOffset-i+this.anchorRect.width,this.horizontalViewportLock&&this.anchorRect.right>this.viewportRect.right&&(this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.right));break;case"insetEnd":this.translateX=this.baseHorizontalOffset,this.horizontalViewportLock&&this.anchorRect.left<this.viewportRect.left&&(this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.left));break;case"end":this.translateX=this.baseHorizontalOffset+this.anchorRect.width,this.horizontalViewportLock&&this.anchorRect.right<this.viewportRect.left&&(this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.left));break;case"center":if(s=(this.anchorRect.width-i)/2,this.translateX=this.baseHorizontalOffset+s,this.horizontalViewportLock){const t=this.anchorRect.left+s,e=this.anchorRect.right-s;t<this.viewportRect.left&&!(e>this.viewportRect.right)?this.translateX=this.translateX-(t-this.viewportRect.left):e>this.viewportRect.right&&!(t<this.viewportRect.left)&&(this.translateX=this.translateX-(e-this.viewportRect.right))}}this.horizontalPosition=t},this.setVerticalPosition=(t,e)=>{if(void 0===t||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect)return;let i=0;switch(this.verticalScaling){case"anchor":case"fill":i=this.verticalViewportLock?this.viewportRect.height:e.height,this.regionHeight=`${i}px`;break;case"content":i=this.regionRect.height,this.regionHeight="unset"}let s=0;switch(t){case"start":this.translateY=this.baseVerticalOffset-i,this.verticalViewportLock&&this.anchorRect.top>this.viewportRect.bottom&&(this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.bottom));break;case"insetStart":this.translateY=this.baseVerticalOffset-i+this.anchorRect.height,this.verticalViewportLock&&this.anchorRect.bottom>this.viewportRect.bottom&&(this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.bottom));break;case"insetEnd":this.translateY=this.baseVerticalOffset,this.verticalViewportLock&&this.anchorRect.top<this.viewportRect.top&&(this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.top));break;case"end":this.translateY=this.baseVerticalOffset+this.anchorRect.height,this.verticalViewportLock&&this.anchorRect.bottom<this.viewportRect.top&&(this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.top));break;case"center":if(s=(this.anchorRect.height-i)/2,this.translateY=this.baseVerticalOffset+s,this.verticalViewportLock){const t=this.anchorRect.top+s,e=this.anchorRect.bottom-s;t<this.viewportRect.top&&!(e>this.viewportRect.bottom)?this.translateY=this.translateY-(t-this.viewportRect.top):e>this.viewportRect.bottom&&!(t<this.viewportRect.top)&&(this.translateY=this.translateY-(e-this.viewportRect.bottom))}}this.verticalPosition=t},this.getPositioningOptions=t=>t?["insetStart","insetEnd"]:["start","end"],this.getAvailableSpace=(t,e,i,s,o,n)=>{const a=e-o,r=n-(e+s);switch(t){case"start":return a;case"insetStart":return a+s;case"insetEnd":return r+s;case"end":return r;case"center":return 2*Math.min(a,r)+s}},this.getNextRegionDimension=(t,e)=>{const i={height:void 0!==this.regionRect?this.regionRect.height:0,width:void 0!==this.regionRect?this.regionRect.width:0};return void 0!==t&&"fill"===this.horizontalScaling?i.width=this.getAvailableSpace(t,void 0!==this.anchorRect?this.anchorRect.left:0,void 0!==this.anchorRect?this.anchorRect.right:0,void 0!==this.anchorRect?this.anchorRect.width:0,void 0!==this.viewportRect?this.viewportRect.left:0,void 0!==this.viewportRect?this.viewportRect.right:0):"anchor"===this.horizontalScaling&&(i.width=void 0!==this.anchorRect?this.anchorRect.width:0),void 0!==e&&"fill"===this.verticalScaling?i.height=this.getAvailableSpace(e,void 0!==this.anchorRect?this.anchorRect.top:0,void 0!==this.anchorRect?this.anchorRect.bottom:0,void 0!==this.anchorRect?this.anchorRect.height:0,void 0!==this.viewportRect?this.viewportRect.top:0,void 0!==this.viewportRect?this.viewportRect.bottom:0):"anchor"===this.verticalScaling&&(i.height=void 0!==this.anchorRect?this.anchorRect.height:0),i},this.startAutoUpdateEventListeners=()=>{window.addEventListener(Vt,this.update,{passive:!0}),window.addEventListener(Nt,this.update,{passive:!0,capture:!0}),null!==this.resizeDetector&&null!==this.viewportElement&&this.resizeDetector.observe(this.viewportElement)},this.stopAutoUpdateEventListeners=()=>{window.removeEventListener(Vt,this.update),window.removeEventListener(Nt,this.update),null!==this.resizeDetector&&null!==this.viewportElement&&this.resizeDetector.unobserve(this.viewportElement)}}anchorChanged(){this.initialLayoutComplete&&(this.anchorElement=this.getAnchor())}viewportChanged(){this.initialLayoutComplete&&(this.viewportElement=this.getViewport())}horizontalPositioningModeChanged(){this.requestReset()}horizontalDefaultPositionChanged(){this.updateForAttributeChange()}horizontalViewportLockChanged(){this.updateForAttributeChange()}horizontalInsetChanged(){this.updateForAttributeChange()}horizontalThresholdChanged(){this.updateForAttributeChange()}horizontalScalingChanged(){this.updateForAttributeChange()}verticalPositioningModeChanged(){this.requestReset()}verticalDefaultPositionChanged(){this.updateForAttributeChange()}verticalViewportLockChanged(){this.updateForAttributeChange()}verticalInsetChanged(){this.updateForAttributeChange()}verticalThresholdChanged(){this.updateForAttributeChange()}verticalScalingChanged(){this.updateForAttributeChange()}fixedPlacementChanged(){this.$fastController.isConnected&&this.initialLayoutComplete&&this.initialize()}autoUpdateModeChanged(t,e){this.$fastController.isConnected&&this.initialLayoutComplete&&("auto"===t&&this.stopAutoUpdateEventListeners(),"auto"===e&&this.startAutoUpdateEventListeners())}anchorElementChanged(){this.requestReset()}viewportElementChanged(){this.$fastController.isConnected&&this.initialLayoutComplete&&this.initialize()}connectedCallback(){super.connectedCallback(),"auto"===this.autoUpdateMode&&this.startAutoUpdateEventListeners(),this.initialize()}disconnectedCallback(){super.disconnectedCallback(),"auto"===this.autoUpdateMode&&this.stopAutoUpdateEventListeners(),this.stopObservers(),this.disconnectResizeDetector()}adoptedCallback(){this.initialize()}disconnectResizeDetector(){null!==this.resizeDetector&&(this.resizeDetector.disconnect(),this.resizeDetector=null)}initializeResizeDetector(){this.disconnectResizeDetector(),this.resizeDetector=new window.ResizeObserver(this.handleResize)}updateForAttributeChange(){this.$fastController.isConnected&&this.initialLayoutComplete&&(this.forceUpdate=!0,this.update())}initialize(){this.initializeResizeDetector(),null===this.anchorElement&&(this.anchorElement=this.getAnchor()),this.requestReset()}requestReset(){this.$fastController.isConnected&&!1===this.pendingReset&&(this.setInitialState(),s.DOM.queueUpdate((()=>this.reset())),this.pendingReset=!0)}setInitialState(){this.initialLayoutComplete=!1,this.regionVisible=!1,this.translateX=0,this.translateY=0,this.baseHorizontalOffset=0,this.baseVerticalOffset=0,this.viewportRect=void 0,this.regionRect=void 0,this.anchorRect=void 0,this.verticalPosition=void 0,this.horizontalPosition=void 0,this.style.opacity="0",this.style.pointerEvents="none",this.forceUpdate=!1,this.style.position=this.fixedPlacement?"fixed":"absolute",this.updatePositionClasses(),this.updateRegionStyle()}}qt.intersectionService=new class{constructor(){this.intersectionDetector=null,this.observedElements=new Map,this.requestPosition=(t,e)=>{var i;null!==this.intersectionDetector&&(this.observedElements.has(t)?null===(i=this.observedElements.get(t))||void 0===i||i.push(e):(this.observedElements.set(t,[e]),this.intersectionDetector.observe(t)))},this.cancelRequestPosition=(t,e)=>{const i=this.observedElements.get(t);if(void 0!==i){const t=i.indexOf(e);-1!==t&&i.splice(t,1)}},this.initializeIntersectionDetector=()=>{s.$global.IntersectionObserver&&(this.intersectionDetector=new IntersectionObserver(this.handleIntersection,{root:null,rootMargin:"0px",threshold:[0,1]}))},this.handleIntersection=t=>{if(null===this.intersectionDetector)return;const e=[],i=[];t.forEach((t=>{var s;null===(s=this.intersectionDetector)||void 0===s||s.unobserve(t.target);const o=this.observedElements.get(t.target);void 0!==o&&(o.forEach((s=>{let o=e.indexOf(s);-1===o&&(o=e.length,e.push(s),i.push([])),i[o].push(t)})),this.observedElements.delete(t.target))})),e.forEach(((t,e)=>{t(i[e])}))},this.initializeIntersectionDetector()}},(0,d.gn)([s.attr],qt.prototype,"anchor",void 0),(0,d.gn)([s.attr],qt.prototype,"viewport",void 0),(0,d.gn)([(0,s.attr)({attribute:"horizontal-positioning-mode"})],qt.prototype,"horizontalPositioningMode",void 0),(0,d.gn)([(0,s.attr)({attribute:"horizontal-default-position"})],qt.prototype,"horizontalDefaultPosition",void 0),(0,d.gn)([(0,s.attr)({attribute:"horizontal-viewport-lock",mode:"boolean"})],qt.prototype,"horizontalViewportLock",void 0),(0,d.gn)([(0,s.attr)({attribute:"horizontal-inset",mode:"boolean"})],qt.prototype,"horizontalInset",void 0),(0,d.gn)([(0,s.attr)({attribute:"horizontal-threshold"})],qt.prototype,"horizontalThreshold",void 0),(0,d.gn)([(0,s.attr)({attribute:"horizontal-scaling"})],qt.prototype,"horizontalScaling",void 0),(0,d.gn)([(0,s.attr)({attribute:"vertical-positioning-mode"})],qt.prototype,"verticalPositioningMode",void 0),(0,d.gn)([(0,s.attr)({attribute:"vertical-default-position"})],qt.prototype,"verticalDefaultPosition",void 0),(0,d.gn)([(0,s.attr)({attribute:"vertical-viewport-lock",mode:"boolean"})],qt.prototype,"verticalViewportLock",void 0),(0,d.gn)([(0,s.attr)({attribute:"vertical-inset",mode:"boolean"})],qt.prototype,"verticalInset",void 0),(0,d.gn)([(0,s.attr)({attribute:"vertical-threshold"})],qt.prototype,"verticalThreshold",void 0),(0,d.gn)([(0,s.attr)({attribute:"vertical-scaling"})],qt.prototype,"verticalScaling",void 0),(0,d.gn)([(0,s.attr)({attribute:"fixed-placement",mode:"boolean"})],qt.prototype,"fixedPlacement",void 0),(0,d.gn)([(0,s.attr)({attribute:"auto-update-mode"})],qt.prototype,"autoUpdateMode",void 0),(0,d.gn)([s.observable],qt.prototype,"anchorElement",void 0),(0,d.gn)([s.observable],qt.prototype,"viewportElement",void 0),(0,d.gn)([s.observable],qt.prototype,"initialLayoutComplete",void 0);const Ut={horizontalDefaultPosition:"center",horizontalPositioningMode:"locktodefault",horizontalInset:!1,horizontalScaling:"anchor"},jt=Object.assign(Object.assign({},Ut),{verticalDefaultPosition:"top",verticalPositioningMode:"locktodefault",verticalInset:!1,verticalScaling:"content"}),_t=Object.assign(Object.assign({},Ut),{verticalDefaultPosition:"bottom",verticalPositioningMode:"locktodefault",verticalInset:!1,verticalScaling:"content"}),Kt=Object.assign(Object.assign({},Ut),{verticalPositioningMode:"dynamic",verticalInset:!1,verticalScaling:"content"}),Wt=Object.assign(Object.assign({},jt),{verticalScaling:"fill"}),Gt=Object.assign(Object.assign({},_t),{verticalScaling:"fill"}),Yt=Object.assign(Object.assign({},Kt),{verticalScaling:"fill"}),Xt=(t,e)=>s.html`
    <div
        class="backplate ${t=>t.shape}"
        part="backplate"
        style="${t=>t.fill?`background-color: var(--avatar-fill-${t.fill});`:void 0}"
    >
        <a
            class="link"
            part="link"
            href="${t=>t.link?t.link:void 0}"
            style="${t=>t.color?`color: var(--avatar-color-${t.color});`:void 0}"
        >
            <slot name="media" part="media">${e.media||""}</slot>
            <slot class="content" part="content"><slot>
        </a>
    </div>
    <slot name="badge" part="badge"></slot>
`;class Qt extends rt{connectedCallback(){super.connectedCallback(),this.shape||(this.shape="circle")}}(0,d.gn)([s.attr],Qt.prototype,"fill",void 0),(0,d.gn)([s.attr],Qt.prototype,"color",void 0),(0,d.gn)([s.attr],Qt.prototype,"link",void 0),(0,d.gn)([s.attr],Qt.prototype,"shape",void 0);const Zt=(t,e)=>s.html`
    <template class="${t=>t.circular?"circular":""}">
        <div class="control" part="control" style="${t=>t.generateBadgeStyle()}">
            <slot></slot>
        </div>
    </template>
`;class Jt extends rt{constructor(){super(...arguments),this.generateBadgeStyle=()=>{if(!this.fill&&!this.color)return;const t=`background-color: var(--badge-fill-${this.fill});`,e=`color: var(--badge-color-${this.color});`;return this.fill&&!this.color?t:this.color&&!this.fill?e:`${e} ${t}`}}}(0,d.gn)([(0,s.attr)({attribute:"fill"})],Jt.prototype,"fill",void 0),(0,d.gn)([(0,s.attr)({attribute:"color"})],Jt.prototype,"color",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Jt.prototype,"circular",void 0);const te=(t,e)=>s.html`
    <div role="listitem" class="listitem" part="listitem">
        ${(0,s.when)((t=>t.href&&t.href.length>0),s.html`
                ${Dt(t,e)}
            `)}
        ${(0,s.when)((t=>!t.href),s.html`
                ${a(t,e)}
                <slot></slot>
                ${n(t,e)}
            `)}
        ${(0,s.when)((t=>t.separator),s.html`
                <span class="separator" part="separator" aria-hidden="true">
                    <slot name="separator">${e.separator||""}</slot>
                </span>
            `)}
    </div>
`;class ee extends St{constructor(){super(...arguments),this.separator=!0}}(0,d.gn)([s.observable],ee.prototype,"separator",void 0),dt(ee,o,At);const ie=(t,e)=>s.html`
    <template role="navigation">
        <div role="list" class="list" part="list">
            <slot
                ${(0,s.slotted)({property:"slottedBreadcrumbItems",filter:(0,s.elements)()})}
            ></slot>
        </div>
    </template>
`;class se extends rt{slottedBreadcrumbItemsChanged(){if(this.$fastController.isConnected){if(void 0===this.slottedBreadcrumbItems||0===this.slottedBreadcrumbItems.length)return;const t=this.slottedBreadcrumbItems[this.slottedBreadcrumbItems.length-1];this.slottedBreadcrumbItems.forEach((e=>{const i=e===t;this.setItemSeparator(e,i),this.setAriaCurrent(e,i)}))}}setItemSeparator(t,e){t instanceof ee&&(t.separator=!e)}findChildWithHref(t){var e,i;return t.childElementCount>0?t.querySelector("a[href]"):(null===(e=t.shadowRoot)||void 0===e?void 0:e.childElementCount)?null===(i=t.shadowRoot)||void 0===i?void 0:i.querySelector("a[href]"):null}setAriaCurrent(t,e){const i=this.findChildWithHref(t);null===i&&t.hasAttribute("href")&&t instanceof ee?e?t.setAttribute("aria-current","page"):t.removeAttribute("aria-current"):null!==i&&(e?i.setAttribute("aria-current","page"):i.removeAttribute("aria-current"))}}(0,d.gn)([s.observable],se.prototype,"slottedBreadcrumbItems",void 0);const oe=(t,e)=>s.html`
    <button
        class="control"
        part="control"
        ?autofocus="${t=>t.autofocus}"
        ?disabled="${t=>t.disabled}"
        form="${t=>t.formId}"
        formaction="${t=>t.formaction}"
        formenctype="${t=>t.formenctype}"
        formmethod="${t=>t.formmethod}"
        formnovalidate="${t=>t.formnovalidate}"
        formtarget="${t=>t.formtarget}"
        name="${t=>t.name}"
        type="${t=>t.type}"
        value="${t=>t.value}"
        aria-atomic="${t=>t.ariaAtomic}"
        aria-busy="${t=>t.ariaBusy}"
        aria-controls="${t=>t.ariaControls}"
        aria-current="${t=>t.ariaCurrent}"
        aria-describedby="${t=>t.ariaDescribedby}"
        aria-details="${t=>t.ariaDetails}"
        aria-disabled="${t=>t.ariaDisabled}"
        aria-errormessage="${t=>t.ariaErrormessage}"
        aria-expanded="${t=>t.ariaExpanded}"
        aria-flowto="${t=>t.ariaFlowto}"
        aria-haspopup="${t=>t.ariaHaspopup}"
        aria-hidden="${t=>t.ariaHidden}"
        aria-invalid="${t=>t.ariaInvalid}"
        aria-keyshortcuts="${t=>t.ariaKeyshortcuts}"
        aria-label="${t=>t.ariaLabel}"
        aria-labelledby="${t=>t.ariaLabelledby}"
        aria-live="${t=>t.ariaLive}"
        aria-owns="${t=>t.ariaOwns}"
        aria-pressed="${t=>t.ariaPressed}"
        aria-relevant="${t=>t.ariaRelevant}"
        aria-roledescription="${t=>t.ariaRoledescription}"
        ${(0,s.ref)("control")}
    >
        ${a(t,e)}
        <span class="content" part="content">
            <slot ${(0,s.slotted)("defaultSlottedContent")}></slot>
        </span>
        ${n(t,e)}
    </button>
`,ne="form-associated-proxy",ae="ElementInternals",re=ae in window&&"setFormValue"in window[ae].prototype,le=new WeakMap;function he(t){const e=class extends t{constructor(...t){super(...t),this.dirtyValue=!1,this.disabled=!1,this.proxyEventsToBlock=["change","click"],this.proxyInitialized=!1,this.required=!1,this.initialValue=this.initialValue||"",this.elementInternals||(this.formResetCallback=this.formResetCallback.bind(this))}static get formAssociated(){return re}get validity(){return this.elementInternals?this.elementInternals.validity:this.proxy.validity}get form(){return this.elementInternals?this.elementInternals.form:this.proxy.form}get validationMessage(){return this.elementInternals?this.elementInternals.validationMessage:this.proxy.validationMessage}get willValidate(){return this.elementInternals?this.elementInternals.willValidate:this.proxy.willValidate}get labels(){if(this.elementInternals)return Object.freeze(Array.from(this.elementInternals.labels));if(this.proxy instanceof HTMLElement&&this.proxy.ownerDocument&&this.id){const t=this.proxy.labels,e=Array.from(this.proxy.getRootNode().querySelectorAll(`[for='${this.id}']`)),i=t?e.concat(Array.from(t)):e;return Object.freeze(i)}return s.emptyArray}valueChanged(t,e){this.dirtyValue=!0,this.proxy instanceof HTMLElement&&(this.proxy.value=this.value),this.currentValue=this.value,this.setFormValue(this.value),this.validate()}currentValueChanged(){this.value=this.currentValue}initialValueChanged(t,e){this.dirtyValue||(this.value=this.initialValue,this.dirtyValue=!1)}disabledChanged(t,e){this.proxy instanceof HTMLElement&&(this.proxy.disabled=this.disabled),s.DOM.queueUpdate((()=>this.classList.toggle("disabled",this.disabled)))}nameChanged(t,e){this.proxy instanceof HTMLElement&&(this.proxy.name=this.name)}requiredChanged(t,e){this.proxy instanceof HTMLElement&&(this.proxy.required=this.required),s.DOM.queueUpdate((()=>this.classList.toggle("required",this.required))),this.validate()}get elementInternals(){if(!re)return null;let t=le.get(this);return t||(t=this.attachInternals(),le.set(this,t)),t}connectedCallback(){super.connectedCallback(),this.addEventListener("keypress",this._keypressHandler),this.value||(this.value=this.initialValue,this.dirtyValue=!1),this.elementInternals||(this.attachProxy(),this.form&&this.form.addEventListener("reset",this.formResetCallback))}disconnectedCallback(){super.disconnectedCallback(),this.proxyEventsToBlock.forEach((t=>this.proxy.removeEventListener(t,this.stopPropagation))),!this.elementInternals&&this.form&&this.form.removeEventListener("reset",this.formResetCallback)}checkValidity(){return this.elementInternals?this.elementInternals.checkValidity():this.proxy.checkValidity()}reportValidity(){return this.elementInternals?this.elementInternals.reportValidity():this.proxy.reportValidity()}setValidity(t,e,i){this.elementInternals?this.elementInternals.setValidity(t,e,i):"string"==typeof e&&this.proxy.setCustomValidity(e)}formDisabledCallback(t){this.disabled=t}formResetCallback(){this.value=this.initialValue,this.dirtyValue=!1}attachProxy(){var t;this.proxyInitialized||(this.proxyInitialized=!0,this.proxy.style.display="none",this.proxyEventsToBlock.forEach((t=>this.proxy.addEventListener(t,this.stopPropagation))),this.proxy.disabled=this.disabled,this.proxy.required=this.required,"string"==typeof this.name&&(this.proxy.name=this.name),"string"==typeof this.value&&(this.proxy.value=this.value),this.proxy.setAttribute("slot",ne),this.proxySlot=document.createElement("slot"),this.proxySlot.setAttribute("name",ne)),null===(t=this.shadowRoot)||void 0===t||t.appendChild(this.proxySlot),this.appendChild(this.proxy)}detachProxy(){var t;this.removeChild(this.proxy),null===(t=this.shadowRoot)||void 0===t||t.removeChild(this.proxySlot)}validate(t){this.proxy instanceof HTMLElement&&this.setValidity(this.proxy.validity,this.proxy.validationMessage,t)}setFormValue(t,e){this.elementInternals&&this.elementInternals.setFormValue(t,e||t)}_keypressHandler(t){if(t.key===ft&&this.form instanceof HTMLFormElement){const t=this.form.querySelector("[type=submit]");null==t||t.click()}}stopPropagation(t){t.stopPropagation()}};return(0,s.attr)({mode:"boolean"})(e.prototype,"disabled"),(0,s.attr)({mode:"fromView",attribute:"value"})(e.prototype,"initialValue"),(0,s.attr)({attribute:"current-value"})(e.prototype,"currentValue"),(0,s.attr)(e.prototype,"name"),(0,s.attr)({mode:"boolean"})(e.prototype,"required"),(0,s.observable)(e.prototype,"value"),e}function de(t){class e extends(he(t)){}class i extends e{constructor(...t){super(t),this.dirtyChecked=!1,this.checkedAttribute=!1,this.checked=!1,this.dirtyChecked=!1}checkedAttributeChanged(){this.defaultChecked=this.checkedAttribute}defaultCheckedChanged(){this.dirtyChecked||(this.checked=this.defaultChecked,this.dirtyChecked=!1)}checkedChanged(t,e){this.dirtyChecked||(this.dirtyChecked=!0),this.currentChecked=this.checked,this.updateForm(),this.proxy instanceof HTMLInputElement&&(this.proxy.checked=this.checked),void 0!==t&&this.$emit("change"),this.validate()}currentCheckedChanged(t,e){this.checked=this.currentChecked}updateForm(){const t=this.checked?this.value:null;this.setFormValue(t,t)}connectedCallback(){super.connectedCallback(),this.updateForm()}formResetCallback(){super.formResetCallback(),this.checked=!!this.checkedAttribute,this.dirtyChecked=!1}}return(0,s.attr)({attribute:"checked",mode:"boolean"})(i.prototype,"checkedAttribute"),(0,s.attr)({attribute:"current-checked",converter:s.booleanConverter})(i.prototype,"currentChecked"),(0,s.observable)(i.prototype,"defaultChecked"),(0,s.observable)(i.prototype,"checked"),i}class ce extends rt{}class ue extends(he(ce)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class pe extends ue{constructor(){super(...arguments),this.handleClick=t=>{var e;this.disabled&&(null===(e=this.defaultSlottedContent)||void 0===e?void 0:e.length)<=1&&t.stopPropagation()},this.handleSubmission=()=>{if(!this.form)return;const t=this.proxy.isConnected;t||this.attachProxy(),"function"==typeof this.form.requestSubmit?this.form.requestSubmit(this.proxy):this.proxy.click(),t||this.detachProxy()},this.handleFormReset=()=>{var t;null===(t=this.form)||void 0===t||t.reset()},this.handleUnsupportedDelegatesFocus=()=>{var t;window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&(null===(t=this.$fastController.definition.shadowOptions)||void 0===t?void 0:t.delegatesFocus)&&(this.focus=()=>{this.control.focus()})}}formactionChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formAction=this.formaction)}formenctypeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formEnctype=this.formenctype)}formmethodChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formMethod=this.formmethod)}formnovalidateChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formNoValidate=this.formnovalidate)}formtargetChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formTarget=this.formtarget)}typeChanged(t,e){this.proxy instanceof HTMLInputElement&&(this.proxy.type=this.type),"submit"===e&&this.addEventListener("click",this.handleSubmission),"submit"===t&&this.removeEventListener("click",this.handleSubmission),"reset"===e&&this.addEventListener("click",this.handleFormReset),"reset"===t&&this.removeEventListener("click",this.handleFormReset)}validate(){super.validate(this.control)}connectedCallback(){var t;super.connectedCallback(),this.proxy.setAttribute("type",this.type),this.handleUnsupportedDelegatesFocus();const e=Array.from(null===(t=this.control)||void 0===t?void 0:t.children);e&&e.forEach((t=>{t.addEventListener("click",this.handleClick)}))}disconnectedCallback(){var t;super.disconnectedCallback();const e=Array.from(null===(t=this.control)||void 0===t?void 0:t.children);e&&e.forEach((t=>{t.removeEventListener("click",this.handleClick)}))}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],pe.prototype,"autofocus",void 0),(0,d.gn)([(0,s.attr)({attribute:"form"})],pe.prototype,"formId",void 0),(0,d.gn)([s.attr],pe.prototype,"formaction",void 0),(0,d.gn)([s.attr],pe.prototype,"formenctype",void 0),(0,d.gn)([s.attr],pe.prototype,"formmethod",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],pe.prototype,"formnovalidate",void 0),(0,d.gn)([s.attr],pe.prototype,"formtarget",void 0),(0,d.gn)([s.attr],pe.prototype,"type",void 0),(0,d.gn)([s.observable],pe.prototype,"defaultSlottedContent",void 0);class me{}(0,d.gn)([(0,s.attr)({attribute:"aria-expanded"})],me.prototype,"ariaExpanded",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-pressed"})],me.prototype,"ariaPressed",void 0),dt(me,Rt),dt(pe,o,me);class ve{constructor(t){if(this.dayFormat="numeric",this.weekdayFormat="long",this.monthFormat="long",this.yearFormat="numeric",this.date=new Date,t)for(const e in t){const i=t[e];"date"===e?this.date=this.getDateObject(i):this[e]=i}}getDateObject(t){if("string"==typeof t){const e=t.split(/[/-]/);return e.length<3?new Date:new Date(parseInt(e[2],10),parseInt(e[0],10)-1,parseInt(e[1],10))}if("day"in t&&"month"in t&&"year"in t){const{day:e,month:i,year:s}=t;return new Date(s,i-1,e)}return t}getDate(t=this.date,e={weekday:this.weekdayFormat,month:this.monthFormat,day:this.dayFormat,year:this.yearFormat},i=this.locale){const s=this.getDateObject(t);if(!s.getTime())return"";const o=Object.assign({timeZone:Intl.DateTimeFormat().resolvedOptions().timeZone},e);return new Intl.DateTimeFormat(i,o).format(s)}getDay(t=this.date.getDate(),e=this.dayFormat,i=this.locale){return this.getDate({month:1,day:t,year:2020},{day:e},i)}getMonth(t=this.date.getMonth()+1,e=this.monthFormat,i=this.locale){return this.getDate({month:t,day:2,year:2020},{month:e},i)}getYear(t=this.date.getFullYear(),e=this.yearFormat,i=this.locale){return this.getDate({month:2,day:2,year:t},{year:e},i)}getWeekday(t=0,e=this.weekdayFormat,i=this.locale){const s=`1-${t+1}-2017`;return this.getDate(s,{weekday:e},i)}getWeekdays(t=this.weekdayFormat,e=this.locale){return Array(7).fill(null).map(((i,s)=>this.getWeekday(s,t,e)))}}class ge extends rt{constructor(){super(...arguments),this.dateFormatter=new ve,this.readonly=!1,this.locale="en-US",this.month=(new Date).getMonth()+1,this.year=(new Date).getFullYear(),this.dayFormat="numeric",this.weekdayFormat="short",this.monthFormat="long",this.yearFormat="numeric",this.minWeeks=0,this.disabledDates="",this.selectedDates="",this.oneDayInMs=864e5}localeChanged(){this.dateFormatter.locale=this.locale}dayFormatChanged(){this.dateFormatter.dayFormat=this.dayFormat}weekdayFormatChanged(){this.dateFormatter.weekdayFormat=this.weekdayFormat}monthFormatChanged(){this.dateFormatter.monthFormat=this.monthFormat}yearFormatChanged(){this.dateFormatter.yearFormat=this.yearFormat}getMonthInfo(t=this.month,e=this.year){const i=t=>new Date(t.getFullYear(),t.getMonth(),1).getDay(),s=t=>{const e=new Date(t.getFullYear(),t.getMonth()+1,1);return new Date(e.getTime()-this.oneDayInMs).getDate()},o=new Date(e,t-1),n=new Date(e,t),a=new Date(e,t-2);return{length:s(o),month:t,start:i(o),year:e,previous:{length:s(a),month:a.getMonth()+1,start:i(a),year:a.getFullYear()},next:{length:s(n),month:n.getMonth()+1,start:i(n),year:n.getFullYear()}}}getDays(t=this.getMonthInfo(),e=this.minWeeks){e=e>10?10:e;const{start:i,length:s,previous:o,next:n}=t,a=[];let r=1-i;for(;r<s+1||a.length<e||a[a.length-1].length%7!=0;){const{month:e,year:i}=r<1?o:r>s?n:t,l=r<1?o.length+r:r>s?r-s:r,h=`${e}-${l}-${i}`,d={day:l,month:e,year:i,disabled:this.dateInString(h,this.disabledDates),selected:this.dateInString(h,this.selectedDates)},c=a[a.length-1];0===a.length||c.length%7==0?a.push([d]):c.push(d),r++}return a}dateInString(t,e){const i=e.split(",").map((t=>t.trim()));return t="string"==typeof t?t:`${t.getMonth()+1}-${t.getDate()}-${t.getFullYear()}`,i.some((e=>e===t))}getDayClassNames(t,e){const{day:i,month:s,year:o,disabled:n,selected:a}=t;return["day",e===`${s}-${i}-${o}`&&"today",this.month!==s&&"inactive",n&&"disabled",a&&"selected"].filter(Boolean).join(" ")}getWeekdayText(){const t=this.dateFormatter.getWeekdays().map((t=>({text:t})));if("long"!==this.weekdayFormat){const e=this.dateFormatter.getWeekdays("long");t.forEach(((t,i)=>{t.abbr=e[i]}))}return t}handleDateSelect(t,e){t.preventDefault,this.$emit("dateselected",e)}handleKeydown(t,e){return t.key===ft&&this.handleDateSelect(t,e),!0}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],ge.prototype,"readonly",void 0),(0,d.gn)([s.attr],ge.prototype,"locale",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],ge.prototype,"month",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],ge.prototype,"year",void 0),(0,d.gn)([(0,s.attr)({attribute:"day-format",mode:"fromView"})],ge.prototype,"dayFormat",void 0),(0,d.gn)([(0,s.attr)({attribute:"weekday-format",mode:"fromView"})],ge.prototype,"weekdayFormat",void 0),(0,d.gn)([(0,s.attr)({attribute:"month-format",mode:"fromView"})],ge.prototype,"monthFormat",void 0),(0,d.gn)([(0,s.attr)({attribute:"year-format",mode:"fromView"})],ge.prototype,"yearFormat",void 0),(0,d.gn)([(0,s.attr)({attribute:"min-weeks",converter:s.nullableNumberConverter})],ge.prototype,"minWeeks",void 0),(0,d.gn)([(0,s.attr)({attribute:"disabled-dates"})],ge.prototype,"disabledDates",void 0),(0,d.gn)([(0,s.attr)({attribute:"selected-dates"})],ge.prototype,"selectedDates",void 0);const be={none:"none",default:"default",sticky:"sticky"},fe={default:"default",columnHeader:"columnheader",rowHeader:"rowheader"},ye={default:"default",header:"header",stickyHeader:"sticky-header"},Ce=s.html`
    <template>
        ${t=>null===t.rowData||null===t.columnDefinition||null===t.columnDefinition.columnDataKey?null:t.rowData[t.columnDefinition.columnDataKey]}
    </template>
`,xe=s.html`
    <template>
        ${t=>null===t.columnDefinition?null:void 0===t.columnDefinition.title?t.columnDefinition.columnDataKey:t.columnDefinition.title}
    </template>
`;class we extends rt{constructor(){super(...arguments),this.cellType=fe.default,this.rowData=null,this.columnDefinition=null,this.isActiveCell=!1,this.customCellView=null,this.updateCellStyle=()=>{this.style.gridColumn=this.gridColumn}}cellTypeChanged(){this.$fastController.isConnected&&this.updateCellView()}gridColumnChanged(){this.$fastController.isConnected&&this.updateCellStyle()}columnDefinitionChanged(t,e){this.$fastController.isConnected&&this.updateCellView()}connectedCallback(){var t;super.connectedCallback(),this.addEventListener(Mt,this.handleFocusin),this.addEventListener(Ht,this.handleFocusout),this.addEventListener(zt,this.handleKeydown),this.style.gridColumn=`${void 0===(null===(t=this.columnDefinition)||void 0===t?void 0:t.gridColumn)?0:this.columnDefinition.gridColumn}`,this.updateCellView(),this.updateCellStyle()}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener(Mt,this.handleFocusin),this.removeEventListener(Ht,this.handleFocusout),this.removeEventListener(zt,this.handleKeydown),this.disconnectCellView()}handleFocusin(t){if(!this.isActiveCell){if(this.isActiveCell=!0,this.cellType===fe.columnHeader){if(null!==this.columnDefinition&&!0!==this.columnDefinition.headerCellInternalFocusQueue&&"function"==typeof this.columnDefinition.headerCellFocusTargetCallback){const t=this.columnDefinition.headerCellFocusTargetCallback(this);null!==t&&t.focus()}}else if(null!==this.columnDefinition&&!0!==this.columnDefinition.cellInternalFocusQueue&&"function"==typeof this.columnDefinition.cellFocusTargetCallback){const t=this.columnDefinition.cellFocusTargetCallback(this);null!==t&&t.focus()}this.$emit("cell-focused",this)}}handleFocusout(t){this===document.activeElement||this.contains(document.activeElement)||(this.isActiveCell=!1)}handleKeydown(t){if(!(t.defaultPrevented||null===this.columnDefinition||this.cellType===fe.default&&!0!==this.columnDefinition.cellInternalFocusQueue||this.cellType===fe.columnHeader&&!0!==this.columnDefinition.headerCellInternalFocusQueue))switch(t.key){case ft:case"F2":if(this.contains(document.activeElement)&&document.activeElement!==this)return;if(this.cellType===fe.columnHeader){if(void 0!==this.columnDefinition.headerCellFocusTargetCallback){const e=this.columnDefinition.headerCellFocusTargetCallback(this);null!==e&&e.focus(),t.preventDefault()}}else if(void 0!==this.columnDefinition.cellFocusTargetCallback){const e=this.columnDefinition.cellFocusTargetCallback(this);null!==e&&e.focus(),t.preventDefault()}break;case yt:this.contains(document.activeElement)&&document.activeElement!==this&&(this.focus(),t.preventDefault())}}updateCellView(){if(this.disconnectCellView(),null!==this.columnDefinition)switch(this.cellType){case fe.columnHeader:void 0!==this.columnDefinition.headerCellTemplate?this.customCellView=this.columnDefinition.headerCellTemplate.render(this,this):this.customCellView=xe.render(this,this);break;case void 0:case fe.rowHeader:case fe.default:void 0!==this.columnDefinition.cellTemplate?this.customCellView=this.columnDefinition.cellTemplate.render(this,this):this.customCellView=Ce.render(this,this)}}disconnectCellView(){null!==this.customCellView&&(this.customCellView.dispose(),this.customCellView=null)}}(0,d.gn)([(0,s.attr)({attribute:"cell-type"})],we.prototype,"cellType",void 0),(0,d.gn)([(0,s.attr)({attribute:"grid-column"})],we.prototype,"gridColumn",void 0),(0,d.gn)([s.observable],we.prototype,"rowData",void 0),(0,d.gn)([s.observable],we.prototype,"columnDefinition",void 0);class $e extends rt{constructor(){super(...arguments),this.rowType=ye.default,this.rowData=null,this.columnDefinitions=null,this.isActiveRow=!1,this.cellsRepeatBehavior=null,this.cellsPlaceholder=null,this.focusColumnIndex=0,this.refocusOnLoad=!1,this.updateRowStyle=()=>{this.style.gridTemplateColumns=this.gridTemplateColumns}}gridTemplateColumnsChanged(){this.$fastController.isConnected&&this.updateRowStyle()}rowTypeChanged(){this.$fastController.isConnected&&this.updateItemTemplate()}rowDataChanged(){null!==this.rowData&&this.isActiveRow&&(this.refocusOnLoad=!0)}cellItemTemplateChanged(){this.updateItemTemplate()}headerCellItemTemplateChanged(){this.updateItemTemplate()}connectedCallback(){super.connectedCallback(),null===this.cellsRepeatBehavior&&(this.cellsPlaceholder=document.createComment(""),this.appendChild(this.cellsPlaceholder),this.updateItemTemplate(),this.cellsRepeatBehavior=new s.RepeatDirective((t=>t.columnDefinitions),(t=>t.activeCellItemTemplate),{positioning:!0}).createBehavior(this.cellsPlaceholder),this.$fastController.addBehaviors([this.cellsRepeatBehavior])),this.addEventListener("cell-focused",this.handleCellFocus),this.addEventListener(Ht,this.handleFocusout),this.addEventListener(zt,this.handleKeydown),this.updateRowStyle(),this.refocusOnLoad&&(this.refocusOnLoad=!1,this.cellElements.length>this.focusColumnIndex&&this.cellElements[this.focusColumnIndex].focus())}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener("cell-focused",this.handleCellFocus),this.removeEventListener(Ht,this.handleFocusout),this.removeEventListener(zt,this.handleKeydown)}handleFocusout(t){this.contains(t.target)||(this.isActiveRow=!1,this.focusColumnIndex=0)}handleCellFocus(t){this.isActiveRow=!0,this.focusColumnIndex=this.cellElements.indexOf(t.target),this.$emit("row-focused",this)}handleKeydown(t){if(t.defaultPrevented)return;let e=0;switch(t.key){case vt:e=Math.max(0,this.focusColumnIndex-1),this.cellElements[e].focus(),t.preventDefault();break;case gt:e=Math.min(this.cellElements.length-1,this.focusColumnIndex+1),this.cellElements[e].focus(),t.preventDefault();break;case Ct:t.ctrlKey||(this.cellElements[0].focus(),t.preventDefault());break;case xt:t.ctrlKey||(this.cellElements[this.cellElements.length-1].focus(),t.preventDefault())}}updateItemTemplate(){this.activeCellItemTemplate=this.rowType===ye.default&&void 0!==this.cellItemTemplate?this.cellItemTemplate:this.rowType===ye.default&&void 0===this.cellItemTemplate?this.defaultCellItemTemplate:void 0!==this.headerCellItemTemplate?this.headerCellItemTemplate:this.defaultHeaderCellItemTemplate}}(0,d.gn)([(0,s.attr)({attribute:"grid-template-columns"})],$e.prototype,"gridTemplateColumns",void 0),(0,d.gn)([(0,s.attr)({attribute:"row-type"})],$e.prototype,"rowType",void 0),(0,d.gn)([s.observable],$e.prototype,"rowData",void 0),(0,d.gn)([s.observable],$e.prototype,"columnDefinitions",void 0),(0,d.gn)([s.observable],$e.prototype,"cellItemTemplate",void 0),(0,d.gn)([s.observable],$e.prototype,"headerCellItemTemplate",void 0),(0,d.gn)([s.observable],$e.prototype,"rowIndex",void 0),(0,d.gn)([s.observable],$e.prototype,"isActiveRow",void 0),(0,d.gn)([s.observable],$e.prototype,"activeCellItemTemplate",void 0),(0,d.gn)([s.observable],$e.prototype,"defaultCellItemTemplate",void 0),(0,d.gn)([s.observable],$e.prototype,"defaultHeaderCellItemTemplate",void 0),(0,d.gn)([s.observable],$e.prototype,"cellElements",void 0);class Ie extends rt{constructor(){super(),this.noTabbing=!1,this.generateHeader=be.default,this.rowsData=[],this.columnDefinitions=null,this.focusRowIndex=0,this.focusColumnIndex=0,this.rowsPlaceholder=null,this.generatedHeader=null,this.isUpdatingFocus=!1,this.pendingFocusUpdate=!1,this.rowindexUpdateQueued=!1,this.columnDefinitionsStale=!0,this.generatedGridTemplateColumns="",this.focusOnCell=(t,e,i)=>{if(0===this.rowElements.length)return this.focusRowIndex=0,void(this.focusColumnIndex=0);const s=Math.max(0,Math.min(this.rowElements.length-1,t)),o=this.rowElements[s].querySelectorAll('[role="cell"], [role="gridcell"], [role="columnheader"], [role="rowheader"]'),n=o[Math.max(0,Math.min(o.length-1,e))];i&&this.scrollHeight!==this.clientHeight&&(s<this.focusRowIndex&&this.scrollTop>0||s>this.focusRowIndex&&this.scrollTop<this.scrollHeight-this.clientHeight)&&n.scrollIntoView({block:"center",inline:"center"}),n.focus()},this.onChildListChange=(t,e)=>{t&&t.length&&(t.forEach((t=>{t.addedNodes.forEach((t=>{1===t.nodeType&&"row"===t.getAttribute("role")&&(t.columnDefinitions=this.columnDefinitions)}))})),this.queueRowIndexUpdate())},this.queueRowIndexUpdate=()=>{this.rowindexUpdateQueued||(this.rowindexUpdateQueued=!0,s.DOM.queueUpdate(this.updateRowIndexes))},this.updateRowIndexes=()=>{let t=this.gridTemplateColumns;if(void 0===t){if(""===this.generatedGridTemplateColumns&&this.rowElements.length>0){const t=this.rowElements[0];this.generatedGridTemplateColumns=new Array(t.cellElements.length).fill("1fr").join(" ")}t=this.generatedGridTemplateColumns}this.rowElements.forEach(((e,i)=>{const s=e;s.rowIndex=i,s.gridTemplateColumns=t,this.columnDefinitionsStale&&(s.columnDefinitions=this.columnDefinitions)})),this.rowindexUpdateQueued=!1,this.columnDefinitionsStale=!1}}static generateTemplateColumns(t){let e="";return t.forEach((t=>{e=`${e}${""===e?"":" "}1fr`})),e}noTabbingChanged(){this.$fastController.isConnected&&(this.noTabbing?this.setAttribute("tabIndex","-1"):this.setAttribute("tabIndex",this.contains(document.activeElement)||this===document.activeElement?"-1":"0"))}generateHeaderChanged(){this.$fastController.isConnected&&this.toggleGeneratedHeader()}gridTemplateColumnsChanged(){this.$fastController.isConnected&&this.updateRowIndexes()}rowsDataChanged(){null===this.columnDefinitions&&this.rowsData.length>0&&(this.columnDefinitions=Ie.generateColumns(this.rowsData[0])),this.$fastController.isConnected&&this.toggleGeneratedHeader()}columnDefinitionsChanged(){null!==this.columnDefinitions?(this.generatedGridTemplateColumns=Ie.generateTemplateColumns(this.columnDefinitions),this.$fastController.isConnected&&(this.columnDefinitionsStale=!0,this.queueRowIndexUpdate())):this.generatedGridTemplateColumns=""}headerCellItemTemplateChanged(){this.$fastController.isConnected&&null!==this.generatedHeader&&(this.generatedHeader.headerCellItemTemplate=this.headerCellItemTemplate)}focusRowIndexChanged(){this.$fastController.isConnected&&this.queueFocusUpdate()}focusColumnIndexChanged(){this.$fastController.isConnected&&this.queueFocusUpdate()}connectedCallback(){super.connectedCallback(),void 0===this.rowItemTemplate&&(this.rowItemTemplate=this.defaultRowItemTemplate),this.rowsPlaceholder=document.createComment(""),this.appendChild(this.rowsPlaceholder),this.toggleGeneratedHeader(),this.rowsRepeatBehavior=new s.RepeatDirective((t=>t.rowsData),(t=>t.rowItemTemplate),{positioning:!0}).createBehavior(this.rowsPlaceholder),this.$fastController.addBehaviors([this.rowsRepeatBehavior]),this.addEventListener("row-focused",this.handleRowFocus),this.addEventListener(Pt,this.handleFocus),this.addEventListener(zt,this.handleKeydown),this.addEventListener(Ht,this.handleFocusOut),this.observer=new MutationObserver(this.onChildListChange),this.observer.observe(this,{childList:!0}),this.noTabbing&&this.setAttribute("tabindex","-1"),s.DOM.queueUpdate(this.queueRowIndexUpdate)}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener("row-focused",this.handleRowFocus),this.removeEventListener(Pt,this.handleFocus),this.removeEventListener(zt,this.handleKeydown),this.removeEventListener(Ht,this.handleFocusOut),this.observer.disconnect(),this.rowsPlaceholder=null,this.generatedHeader=null}handleRowFocus(t){this.isUpdatingFocus=!0;const e=t.target;this.focusRowIndex=this.rowElements.indexOf(e),this.focusColumnIndex=e.focusColumnIndex,this.setAttribute("tabIndex","-1"),this.isUpdatingFocus=!1}handleFocus(t){this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,!0)}handleFocusOut(t){null!==t.relatedTarget&&this.contains(t.relatedTarget)||this.setAttribute("tabIndex",this.noTabbing?"-1":"0")}handleKeydown(t){if(t.defaultPrevented)return;let e;const i=this.rowElements.length-1,s=this.offsetHeight+this.scrollTop,o=this.rowElements[i];switch(t.key){case bt:t.preventDefault(),this.focusOnCell(this.focusRowIndex-1,this.focusColumnIndex,!0);break;case mt:t.preventDefault(),this.focusOnCell(this.focusRowIndex+1,this.focusColumnIndex,!0);break;case"PageUp":if(t.preventDefault(),0===this.rowElements.length){this.focusOnCell(0,0,!1);break}if(0===this.focusRowIndex)return void this.focusOnCell(0,this.focusColumnIndex,!1);for(e=this.focusRowIndex-1;e>=0;e--){const t=this.rowElements[e];if(t.offsetTop<this.scrollTop){this.scrollTop=t.offsetTop+t.clientHeight-this.clientHeight;break}}this.focusOnCell(e,this.focusColumnIndex,!1);break;case"PageDown":if(t.preventDefault(),0===this.rowElements.length){this.focusOnCell(0,0,!1);break}if(this.focusRowIndex>=i||o.offsetTop+o.offsetHeight<=s)return void this.focusOnCell(i,this.focusColumnIndex,!1);for(e=this.focusRowIndex+1;e<=i;e++){const t=this.rowElements[e];if(t.offsetTop+t.offsetHeight>s){let e=0;this.generateHeader===be.sticky&&null!==this.generatedHeader&&(e=this.generatedHeader.clientHeight),this.scrollTop=t.offsetTop-e;break}}this.focusOnCell(e,this.focusColumnIndex,!1);break;case Ct:t.ctrlKey&&(t.preventDefault(),this.focusOnCell(0,0,!0));break;case xt:t.ctrlKey&&null!==this.columnDefinitions&&(t.preventDefault(),this.focusOnCell(this.rowElements.length-1,this.columnDefinitions.length-1,!0))}}queueFocusUpdate(){this.isUpdatingFocus&&(this.contains(document.activeElement)||this===document.activeElement)||!1===this.pendingFocusUpdate&&(this.pendingFocusUpdate=!0,s.DOM.queueUpdate((()=>this.updateFocus())))}updateFocus(){this.pendingFocusUpdate=!1,this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,!0)}toggleGeneratedHeader(){if(null!==this.generatedHeader&&(this.removeChild(this.generatedHeader),this.generatedHeader=null),this.generateHeader!==be.none&&this.rowsData.length>0){const t=document.createElement(this.rowElementTag);return this.generatedHeader=t,this.generatedHeader.columnDefinitions=this.columnDefinitions,this.generatedHeader.gridTemplateColumns=this.gridTemplateColumns,this.generatedHeader.rowType=this.generateHeader===be.sticky?ye.stickyHeader:ye.header,void(null===this.firstChild&&null===this.rowsPlaceholder||this.insertBefore(t,null!==this.firstChild?this.firstChild:this.rowsPlaceholder))}}}Ie.generateColumns=t=>Object.getOwnPropertyNames(t).map(((t,e)=>({columnDataKey:t,gridColumn:`${e}`}))),(0,d.gn)([(0,s.attr)({attribute:"no-tabbing",mode:"boolean"})],Ie.prototype,"noTabbing",void 0),(0,d.gn)([(0,s.attr)({attribute:"generate-header"})],Ie.prototype,"generateHeader",void 0),(0,d.gn)([(0,s.attr)({attribute:"grid-template-columns"})],Ie.prototype,"gridTemplateColumns",void 0),(0,d.gn)([s.observable],Ie.prototype,"rowsData",void 0),(0,d.gn)([s.observable],Ie.prototype,"columnDefinitions",void 0),(0,d.gn)([s.observable],Ie.prototype,"rowItemTemplate",void 0),(0,d.gn)([s.observable],Ie.prototype,"cellItemTemplate",void 0),(0,d.gn)([s.observable],Ie.prototype,"headerCellItemTemplate",void 0),(0,d.gn)([s.observable],Ie.prototype,"focusRowIndex",void 0),(0,d.gn)([s.observable],Ie.prototype,"focusColumnIndex",void 0),(0,d.gn)([s.observable],Ie.prototype,"defaultRowItemTemplate",void 0),(0,d.gn)([s.observable],Ie.prototype,"rowElementTag",void 0),(0,d.gn)([s.observable],Ie.prototype,"rowElements",void 0);const ke=s.html`
    <div
        class="title"
        part="title"
        aria-label="${t=>t.dateFormatter.getDate(`${t.month}-2-${t.year}`,{month:"long",year:"numeric"})}"
    >
        <span part="month">
            ${t=>t.dateFormatter.getMonth(t.month)}
        </span>
        <span part="year">${t=>t.dateFormatter.getYear(t.year)}</span>
    </div>
`,Ee=t=>{const e=t.tagFor(we);return s.html`
        <${e}
            class="week-day"
            part="week-day"
            tabindex="-1"
            grid-column="${(t,e)=>e.index+1}"
            abbr="${t=>t.abbr}"
        >
            ${t=>t.text}
        </${e}>
    `},Te=(t,e)=>{const i=t.tagFor(we);return s.html`
        <${i}
            class="${(t,i)=>i.parentContext.parent.getDayClassNames(t,e)}"
            part="day"
            tabindex="-1"
            role="gridcell"
            grid-column="${(t,e)=>e.index+1}"
            @click="${(t,e)=>e.parentContext.parent.handleDateSelect(e.event,t)}"
            @keydown="${(t,e)=>e.parentContext.parent.handleKeydown(e.event,t)}"
            aria-label="${(t,e)=>e.parentContext.parent.dateFormatter.getDate(`${t.month}-${t.day}-${t.year}`,{month:"long",day:"numeric"})}"
        >
            <div
                class="date"
                part="${t=>e===`${t.month}-${t.day}-${t.year}`?"today":"date"}"
            >
                ${(t,e)=>e.parentContext.parent.dateFormatter.getDay(t.day)}
            </div>
            <slot name="${t=>t.month}-${t=>t.day}-${t=>t.year}"></slot>
        </${i}>
    `},Oe=(t,e)=>{const i=t.tagFor($e);return s.html`
        <${i}
            class="week"
            part="week"
            role="row"
            role-type="default"
            grid-template-columns="1fr 1fr 1fr 1fr 1fr 1fr 1fr"
        >
        ${(0,s.repeat)((t=>t),Te(t,e),{positioning:!0})}
        </${i}>
    `},De=(t,e)=>{const i=t.tagFor(Ie),o=t.tagFor($e);return s.html`
    <${i} class="days interact" part="days" generate-header="none">
        <${o}
            class="week-days"
            part="week-days"
            role="row"
            row-type="header"
            grid-template-columns="1fr 1fr 1fr 1fr 1fr 1fr 1fr"
        >
            ${(0,s.repeat)((t=>t.getWeekdayText()),Ee(t),{positioning:!0})}
        </${o}>
        ${(0,s.repeat)((t=>t.getDays()),Oe(t,e))}
    </${i}>
`},Re=t=>s.html`
        <div class="days" part="days">
            <div class="week-days" part="week-days">
                ${(0,s.repeat)((t=>t.getWeekdayText()),s.html`
                        <div class="week-day" part="week-day" abbr="${t=>t.abbr}">
                            ${t=>t.text}
                        </div>
                    `)}
            </div>
            ${(0,s.repeat)((t=>t.getDays()),s.html`
                    <div class="week">
                        ${(0,s.repeat)((t=>t),s.html`
                                <div
                                    class="${(e,i)=>i.parentContext.parent.getDayClassNames(e,t)}"
                                    part="day"
                                    aria-label="${(t,e)=>e.parentContext.parent.dateFormatter.getDate(`${t.month}-${t.day}-${t.year}`,{month:"long",day:"numeric"})}"
                                >
                                    <div
                                        class="date"
                                        part="${e=>t===`${e.month}-${e.day}-${e.year}`?"today":"date"}"
                                    >
                                        ${(t,e)=>e.parentContext.parent.dateFormatter.getDay(t.day)}
                                    </div>
                                    <slot
                                        name="${t=>t.month}-${t=>t.day}-${t=>t.year}"
                                    ></slot>
                                </div>
                            `)}
                    </div>
                `)}
        </div>
    `,Se=(t,e)=>{var i;const o=new Date,n=`${o.getMonth()+1}-${o.getDate()}-${o.getFullYear()}`;return s.html`
        <template>
            ${l}
            ${e.title instanceof Function?e.title(t,e):null!==(i=e.title)&&void 0!==i?i:""}
            <slot></slot>
            ${(0,s.when)((t=>t.readonly),Re(n),De(t,n))}
            ${r}
        </template>
    `},Ae=(t,e)=>s.html`
    <slot></slot>
`;class Fe extends rt{}const Le=(t,e)=>s.html`
    <template
        role="checkbox"
        aria-checked="${t=>t.checked}"
        aria-required="${t=>t.required}"
        aria-disabled="${t=>t.disabled}"
        aria-readonly="${t=>t.readOnly}"
        tabindex="${t=>t.disabled?null:0}"
        @keypress="${(t,e)=>t.keypressHandler(e.event)}"
        @click="${(t,e)=>t.clickHandler(e.event)}"
        class="${t=>t.readOnly?"readonly":""} ${t=>t.checked?"checked":""} ${t=>t.indeterminate?"indeterminate":""}"
    >
        <div part="control" class="control">
            <slot name="checked-indicator">
                ${e.checkedIndicator||""}
            </slot>
            <slot name="indeterminate-indicator">
                ${e.indeterminateIndicator||""}
            </slot>
        </div>
        <label
            part="label"
            class="${t=>t.defaultSlottedNodes&&t.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
    </template>
`;class Pe extends rt{}class Me extends(de(Pe)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class He extends Me{constructor(){super(),this.initialValue="on",this.indeterminate=!1,this.keypressHandler=t=>{this.readOnly||t.key!==wt||(this.indeterminate&&(this.indeterminate=!1),this.checked=!this.checked)},this.clickHandler=t=>{this.disabled||this.readOnly||(this.indeterminate&&(this.indeterminate=!1),this.checked=!this.checked)},this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],He.prototype,"readOnly",void 0),(0,d.gn)([s.observable],He.prototype,"defaultSlottedNodes",void 0),(0,d.gn)([s.observable],He.prototype,"indeterminate",void 0);let ze,Ve=0;function Ne(t=""){return`${t}${Ve++}`}function Be(...t){return t.every((t=>t instanceof HTMLElement))}function qe(t){return Be(t)&&("option"===t.getAttribute("role")||t instanceof HTMLOptionElement)}class Ue extends rt{constructor(t,e,i,s){super(),this.defaultSelected=!1,this.dirtySelected=!1,this.selected=this.defaultSelected,this.dirtyValue=!1,t&&(this.textContent=t),e&&(this.initialValue=e),i&&(this.defaultSelected=i),s&&(this.selected=s),this.proxy=new Option(`${this.textContent}`,this.initialValue,this.defaultSelected,this.selected),this.proxy.disabled=this.disabled}checkedChanged(t,e){this.ariaChecked="boolean"!=typeof e?null:e?"true":"false"}contentChanged(t,e){this.proxy instanceof HTMLOptionElement&&(this.proxy.textContent=this.textContent),this.$emit("contentchange",null,{bubbles:!0})}defaultSelectedChanged(){this.dirtySelected||(this.selected=this.defaultSelected,this.proxy instanceof HTMLOptionElement&&(this.proxy.selected=this.defaultSelected))}disabledChanged(t,e){this.ariaDisabled=this.disabled?"true":"false",this.proxy instanceof HTMLOptionElement&&(this.proxy.disabled=this.disabled)}selectedAttributeChanged(){this.defaultSelected=this.selectedAttribute,this.proxy instanceof HTMLOptionElement&&(this.proxy.defaultSelected=this.defaultSelected)}selectedChanged(){this.ariaSelected=this.selected?"true":"false",this.dirtySelected||(this.dirtySelected=!0),this.proxy instanceof HTMLOptionElement&&(this.proxy.selected=this.selected)}initialValueChanged(t,e){this.dirtyValue||(this.value=this.initialValue,this.dirtyValue=!1)}get label(){var t;return null!==(t=this.value)&&void 0!==t?t:this.text}get text(){var t,e;return null!==(e=null===(t=this.textContent)||void 0===t?void 0:t.replace(/\s+/g," ").trim())&&void 0!==e?e:""}set value(t){const e=`${null!=t?t:""}`;this._value=e,this.dirtyValue=!0,this.proxy instanceof HTMLOptionElement&&(this.proxy.value=e),s.Observable.notify(this,"value")}get value(){var t;return s.Observable.track(this,"value"),null!==(t=this._value)&&void 0!==t?t:this.text}get form(){return this.proxy?this.proxy.form:null}}(0,d.gn)([s.observable],Ue.prototype,"checked",void 0),(0,d.gn)([s.observable],Ue.prototype,"content",void 0),(0,d.gn)([s.observable],Ue.prototype,"defaultSelected",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Ue.prototype,"disabled",void 0),(0,d.gn)([(0,s.attr)({attribute:"selected",mode:"boolean"})],Ue.prototype,"selectedAttribute",void 0),(0,d.gn)([s.observable],Ue.prototype,"selected",void 0),(0,d.gn)([(0,s.attr)({attribute:"value",mode:"fromView"})],Ue.prototype,"initialValue",void 0);class je{}(0,d.gn)([s.observable],je.prototype,"ariaChecked",void 0),(0,d.gn)([s.observable],je.prototype,"ariaPosInSet",void 0),(0,d.gn)([s.observable],je.prototype,"ariaSelected",void 0),(0,d.gn)([s.observable],je.prototype,"ariaSetSize",void 0),dt(je,Rt),dt(Ue,o,je);class _e extends rt{constructor(){super(...arguments),this._options=[],this.selectedIndex=-1,this.selectedOptions=[],this.shouldSkipFocus=!1,this.typeaheadBuffer="",this.typeaheadExpired=!0,this.typeaheadTimeout=-1}get firstSelectedOption(){var t;return null!==(t=this.selectedOptions[0])&&void 0!==t?t:null}get hasSelectableOptions(){return this.options.length>0&&!this.options.every((t=>t.disabled))}get length(){var t,e;return null!==(e=null===(t=this.options)||void 0===t?void 0:t.length)&&void 0!==e?e:0}get options(){return s.Observable.track(this,"options"),this._options}set options(t){this._options=t,s.Observable.notify(this,"options")}get typeAheadExpired(){return this.typeaheadExpired}set typeAheadExpired(t){this.typeaheadExpired=t}clickHandler(t){const e=t.target.closest("option,[role=option]");if(e&&!e.disabled)return this.selectedIndex=this.options.indexOf(e),!0}focusAndScrollOptionIntoView(t=this.firstSelectedOption){this.contains(document.activeElement)&&null!==t&&(t.focus(),requestAnimationFrame((()=>{t.scrollIntoView({block:"nearest"})})))}focusinHandler(t){this.shouldSkipFocus||t.target!==t.currentTarget||(this.setSelectedOptions(),this.focusAndScrollOptionIntoView()),this.shouldSkipFocus=!1}getTypeaheadMatches(){const t=this.typeaheadBuffer.replace(/[.*+\-?^${}()|[\]\\]/g,"\\$&"),e=new RegExp(`^${t}`,"gi");return this.options.filter((t=>t.text.trim().match(e)))}getSelectableIndex(t=this.selectedIndex,e){const i=t>e?-1:t<e?1:0,s=t+i;let o=null;switch(i){case-1:o=this.options.reduceRight(((t,e,i)=>!t&&!e.disabled&&i<s?e:t),o);break;case 1:o=this.options.reduce(((t,e,i)=>!t&&!e.disabled&&i>s?e:t),o)}return this.options.indexOf(o)}handleChange(t,e){"selected"===e&&(_e.slottedOptionFilter(t)&&(this.selectedIndex=this.options.indexOf(t)),this.setSelectedOptions())}handleTypeAhead(t){this.typeaheadTimeout&&window.clearTimeout(this.typeaheadTimeout),this.typeaheadTimeout=window.setTimeout((()=>this.typeaheadExpired=!0),_e.TYPE_AHEAD_TIMEOUT_MS),t.length>1||(this.typeaheadBuffer=`${this.typeaheadExpired?"":this.typeaheadBuffer}${t}`)}keydownHandler(t){if(this.disabled)return!0;this.shouldSkipFocus=!1;const e=t.key;switch(e){case Ct:t.shiftKey||(t.preventDefault(),this.selectFirstOption());break;case mt:t.shiftKey||(t.preventDefault(),this.selectNextOption());break;case bt:t.shiftKey||(t.preventDefault(),this.selectPreviousOption());break;case xt:t.preventDefault(),this.selectLastOption();break;case $t:return this.focusAndScrollOptionIntoView(),!0;case ft:case yt:return!0;case wt:if(this.typeaheadExpired)return!0;default:return 1===e.length&&this.handleTypeAhead(`${e}`),!0}}mousedownHandler(t){return this.shouldSkipFocus=!this.contains(document.activeElement),!0}multipleChanged(t,e){this.ariaMultiSelectable=e?"true":null}selectedIndexChanged(t,e){var i;if(this.hasSelectableOptions){if((null===(i=this.options[this.selectedIndex])||void 0===i?void 0:i.disabled)&&"number"==typeof t){const i=this.getSelectableIndex(t,e),s=i>-1?i:t;return this.selectedIndex=s,void(e===s&&this.selectedIndexChanged(e,s))}this.setSelectedOptions()}else this.selectedIndex=-1}selectedOptionsChanged(t,e){var i;const o=e.filter(_e.slottedOptionFilter);null===(i=this.options)||void 0===i||i.forEach((t=>{const e=s.Observable.getNotifier(t);e.unsubscribe(this,"selected"),t.selected=o.includes(t),e.subscribe(this,"selected")}))}selectFirstOption(){var t,e;this.disabled||(this.selectedIndex=null!==(e=null===(t=this.options)||void 0===t?void 0:t.findIndex((t=>!t.disabled)))&&void 0!==e?e:-1)}selectLastOption(){this.disabled||(this.selectedIndex=function(t,e){let i=t.length;for(;i--;)if(!t[i].disabled)return i;return-1}(this.options))}selectNextOption(){!this.disabled&&this.selectedIndex<this.options.length-1&&(this.selectedIndex+=1)}selectPreviousOption(){!this.disabled&&this.selectedIndex>0&&(this.selectedIndex=this.selectedIndex-1)}setDefaultSelectedOption(){var t,e;this.selectedIndex=null!==(e=null===(t=this.options)||void 0===t?void 0:t.findIndex((t=>t.defaultSelected)))&&void 0!==e?e:-1}setSelectedOptions(){var t,e,i;(null===(t=this.options)||void 0===t?void 0:t.length)&&(this.selectedOptions=[this.options[this.selectedIndex]],this.ariaActiveDescendant=null!==(i=null===(e=this.firstSelectedOption)||void 0===e?void 0:e.id)&&void 0!==i?i:"",this.focusAndScrollOptionIntoView())}slottedOptionsChanged(t,e){this.options=e.reduce(((t,e)=>(qe(e)&&t.push(e),t)),[]);const i=`${this.options.length}`;this.options.forEach(((t,e)=>{t.id||(t.id=Ne("option-")),t.ariaPosInSet=`${e+1}`,t.ariaSetSize=i})),this.$fastController.isConnected&&(this.setSelectedOptions(),this.setDefaultSelectedOption())}typeaheadBufferChanged(t,e){if(this.$fastController.isConnected){const t=this.getTypeaheadMatches();if(t.length){const e=this.options.indexOf(t[0]);e>-1&&(this.selectedIndex=e)}this.typeaheadExpired=!1}}}_e.slottedOptionFilter=t=>qe(t)&&!t.hidden,_e.TYPE_AHEAD_TIMEOUT_MS=1e3,(0,d.gn)([(0,s.attr)({mode:"boolean"})],_e.prototype,"disabled",void 0),(0,d.gn)([s.observable],_e.prototype,"selectedIndex",void 0),(0,d.gn)([s.observable],_e.prototype,"selectedOptions",void 0),(0,d.gn)([s.observable],_e.prototype,"slottedOptions",void 0),(0,d.gn)([s.observable],_e.prototype,"typeaheadBuffer",void 0);class Ke{}(0,d.gn)([s.observable],Ke.prototype,"ariaActiveDescendant",void 0),(0,d.gn)([s.observable],Ke.prototype,"ariaDisabled",void 0),(0,d.gn)([s.observable],Ke.prototype,"ariaExpanded",void 0),(0,d.gn)([s.observable],Ke.prototype,"ariaMultiSelectable",void 0),dt(Ke,Rt),dt(_e,Ke);const We={above:"above",below:"below"};class Ge extends _e{}class Ye extends(he(Ge)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}const Xe={inline:"inline",list:"list",both:"both",none:"none"};class Qe extends Ye{constructor(){super(...arguments),this._value="",this.filteredOptions=[],this.filter="",this.forcedPosition=!1,this.listboxId=Ne("listbox-"),this.maxHeight=0,this.open=!1}formResetCallback(){super.formResetCallback(),this.setDefaultSelectedOption(),this.updateValue()}validate(){super.validate(this.control)}get isAutocompleteInline(){return this.autocomplete===Xe.inline||this.isAutocompleteBoth}get isAutocompleteList(){return this.autocomplete===Xe.list||this.isAutocompleteBoth}get isAutocompleteBoth(){return this.autocomplete===Xe.both}openChanged(){if(this.open)return this.ariaControls=this.listboxId,this.ariaExpanded="true",this.setPositioning(),this.focusAndScrollOptionIntoView(),void s.DOM.queueUpdate((()=>this.focus()));this.ariaControls="",this.ariaExpanded="false"}get options(){return s.Observable.track(this,"options"),this.filteredOptions.length?this.filteredOptions:this._options}set options(t){this._options=t,s.Observable.notify(this,"options")}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}positionChanged(t,e){this.positionAttribute=e,this.setPositioning()}get value(){return s.Observable.track(this,"value"),this._value}set value(t){var e,i,o;const n=`${this._value}`;if(this.$fastController.isConnected&&this.options){const s=this.options.findIndex((e=>e.text.toLowerCase()===t.toLowerCase())),n=null===(e=this.options[this.selectedIndex])||void 0===e?void 0:e.text,a=null===(i=this.options[s])||void 0===i?void 0:i.text;this.selectedIndex=n!==a?s:this.selectedIndex,t=(null===(o=this.firstSelectedOption)||void 0===o?void 0:o.text)||t}n!==t&&(this._value=t,super.valueChanged(n,t),s.Observable.notify(this,"value"))}clickHandler(t){if(!this.disabled){if(this.open){const e=t.target.closest("option,[role=option]");if(!e||e.disabled)return;this.selectedOptions=[e],this.control.value=e.text,this.clearSelectionRange(),this.updateValue(!0)}return this.open=!this.open,this.open&&this.control.focus(),!0}}connectedCallback(){super.connectedCallback(),this.forcedPosition=!!this.positionAttribute,this.value&&(this.initialValue=this.value)}disabledChanged(t,e){super.disabledChanged&&super.disabledChanged(t,e),this.ariaDisabled=this.disabled?"true":"false"}filterOptions(){this.autocomplete&&this.autocomplete!==Xe.none||(this.filter="");const t=this.filter.toLowerCase();this.filteredOptions=this._options.filter((t=>t.text.toLowerCase().startsWith(this.filter.toLowerCase()))),this.isAutocompleteList&&(this.filteredOptions.length||t||(this.filteredOptions=this._options),this._options.forEach((t=>{t.hidden=!this.filteredOptions.includes(t)})))}focusAndScrollOptionIntoView(){this.contains(document.activeElement)&&(this.control.focus(),this.firstSelectedOption&&requestAnimationFrame((()=>{var t;null===(t=this.firstSelectedOption)||void 0===t||t.scrollIntoView({block:"nearest"})})))}focusoutHandler(t){if(this.syncValue(),!this.open)return!0;const e=t.relatedTarget;this.isSameNode(e)?this.focus():this.options&&this.options.includes(e)||(this.open=!1)}inputHandler(t){if(this.filter=this.control.value,this.filterOptions(),this.isAutocompleteInline||(this.selectedIndex=this.options.map((t=>t.text)).indexOf(this.control.value)),t.inputType.includes("deleteContent")||!this.filter.length)return!0;this.isAutocompleteList&&!this.open&&(this.open=!0),this.isAutocompleteInline&&(this.filteredOptions.length?(this.selectedOptions=[this.filteredOptions[0]],this.selectedIndex=this.options.indexOf(this.firstSelectedOption),this.setInlineSelection()):this.selectedIndex=-1)}keydownHandler(t){const e=t.key;if(t.ctrlKey||t.shiftKey)return!0;switch(e){case"Enter":this.syncValue(),this.isAutocompleteInline&&(this.filter=this.value),this.open=!1,this.clearSelectionRange();break;case"Escape":if(this.isAutocompleteInline||(this.selectedIndex=-1),this.open){this.open=!1;break}this.value="",this.control.value="",this.filter="",this.filterOptions();break;case"Tab":if(this.setInputToSelection(),!this.open)return!0;t.preventDefault(),this.open=!1;break;case"ArrowUp":case"ArrowDown":if(this.filterOptions(),!this.open){this.open=!0;break}this.filteredOptions.length>0&&super.keydownHandler(t),this.isAutocompleteInline&&this.setInlineSelection();break;default:return!0}}keyupHandler(t){switch(t.key){case"ArrowLeft":case"ArrowRight":case"Backspace":case"Delete":case"Home":case"End":this.filter=this.control.value,this.selectedIndex=-1,this.filterOptions()}}selectedIndexChanged(t,e){if(this.$fastController.isConnected){if((e=kt(-1,this.options.length-1,e))!==this.selectedIndex)return void(this.selectedIndex=e);super.selectedIndexChanged(t,e)}}selectPreviousOption(){!this.disabled&&this.selectedIndex>=0&&(this.selectedIndex=this.selectedIndex-1)}setDefaultSelectedOption(){if(this.$fastController.isConnected&&this.options){const t=this.options.findIndex((t=>null!==t.getAttribute("selected")||t.selected));this.selectedIndex=t,!this.dirtyValue&&this.firstSelectedOption&&(this.value=this.firstSelectedOption.text),this.setSelectedOptions()}}setInputToSelection(){this.firstSelectedOption&&(this.control.value=this.firstSelectedOption.text,this.control.focus())}setInlineSelection(){this.firstSelectedOption&&(this.setInputToSelection(),this.control.setSelectionRange(this.filter.length,this.control.value.length,"backward"))}syncValue(){var t;const e=this.selectedIndex>-1?null===(t=this.firstSelectedOption)||void 0===t?void 0:t.text:this.control.value;this.updateValue(this.value!==e)}setPositioning(){const t=this.getBoundingClientRect(),e=window.innerHeight-t.bottom;this.position=this.forcedPosition?this.positionAttribute:t.top>e?We.above:We.below,this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position,this.maxHeight=this.position===We.above?~~t.top:~~e}selectedOptionsChanged(t,e){this.$fastController.isConnected&&this._options.forEach((t=>{t.selected=e.includes(t)}))}slottedOptionsChanged(t,e){super.slottedOptionsChanged(t,e),this.updateValue()}updateValue(t){var e;this.$fastController.isConnected&&(this.value=(null===(e=this.firstSelectedOption)||void 0===e?void 0:e.text)||this.control.value,this.control.value=this.value),t&&this.$emit("change")}clearSelectionRange(){const t=this.control.value.length;this.control.setSelectionRange(t,t)}}(0,d.gn)([(0,s.attr)({attribute:"autocomplete",mode:"fromView"})],Qe.prototype,"autocomplete",void 0),(0,d.gn)([s.observable],Qe.prototype,"maxHeight",void 0),(0,d.gn)([(0,s.attr)({attribute:"open",mode:"boolean"})],Qe.prototype,"open",void 0),(0,d.gn)([s.attr],Qe.prototype,"placeholder",void 0),(0,d.gn)([(0,s.attr)({attribute:"position"})],Qe.prototype,"positionAttribute",void 0),(0,d.gn)([s.observable],Qe.prototype,"position",void 0);class Ze{}(0,d.gn)([s.observable],Ze.prototype,"ariaAutoComplete",void 0),(0,d.gn)([s.observable],Ze.prototype,"ariaControls",void 0),dt(Ze,Ke),dt(Qe,o,Ze);const Je=(t,e)=>s.html`
    <template
        aria-disabled="${t=>t.ariaDisabled}"
        autocomplete="${t=>t.autocomplete}"
        class="${t=>t.open?"open":""} ${t=>t.disabled?"disabled":""} ${t=>t.position}"
        ?open="${t=>t.open}"
        tabindex="${t=>t.disabled?null:"0"}"
        @click="${(t,e)=>t.clickHandler(e.event)}"
        @focusout="${(t,e)=>t.focusoutHandler(e.event)}"
        @keydown="${(t,e)=>t.keydownHandler(e.event)}"
    >
        <div class="control" part="control">
            ${a(t,e)}
            <slot name="control">
                <input
                    aria-activedescendant="${t=>t.open?t.ariaActiveDescendant:null}"
                    aria-autocomplete="${t=>t.ariaAutoComplete}"
                    aria-controls="${t=>t.ariaControls}"
                    aria-disabled="${t=>t.ariaDisabled}"
                    aria-expanded="${t=>t.ariaExpanded}"
                    aria-haspopup="listbox"
                    class="selected-value"
                    part="selected-value"
                    placeholder="${t=>t.placeholder}"
                    role="combobox"
                    type="text"
                    ?disabled="${t=>t.disabled}"
                    :value="${t=>t.value}"
                    @input="${(t,e)=>t.inputHandler(e.event)}"
                    @keyup="${(t,e)=>t.keyupHandler(e.event)}"
                    ${(0,s.ref)("control")}
                />
                <div class="indicator" part="indicator" aria-hidden="true">
                    <slot name="indicator">
                        ${e.indicator||""}
                    </slot>
                </div>
            </slot>
            ${n(t,e)}
        </div>
        <div
            class="listbox"
            id="${t=>t.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${t=>t.disabled}"
            ?hidden="${t=>!t.open}"
            ${(0,s.ref)("listbox")}
        >
            <slot
                ${(0,s.slotted)({filter:_e.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`,ti=(t,e)=>{const i=function(t){const e=t.tagFor($e);return s.html`
    <${e}
        :rowData="${t=>t}"
        :cellItemTemplate="${(t,e)=>e.parent.cellItemTemplate}"
        :headerCellItemTemplate="${(t,e)=>e.parent.headerCellItemTemplate}"
    ></${e}>
`}(t),o=t.tagFor($e);return s.html`
        <template
            role="grid"
            tabindex="0"
            :rowElementTag="${()=>o}"
            :defaultRowItemTemplate="${i}"
            ${(0,s.children)({property:"rowElements",filter:(0,s.elements)("[role=row]")})}
        >
            <slot></slot>
        </template>
    `},ei=(t,e)=>{const i=function(t){const e=t.tagFor(we);return s.html`
    <${e}
        cell-type="${t=>t.isRowHeader?"rowheader":void 0}"
        grid-column="${(t,e)=>e.index+1}"
        :rowData="${(t,e)=>e.parent.rowData}"
        :columnDefinition="${t=>t}"
    ></${e}>
`}(t),o=function(t){const e=t.tagFor(we);return s.html`
    <${e}
        cell-type="columnheader"
        grid-column="${(t,e)=>e.index+1}"
        :columnDefinition="${t=>t}"
    ></${e}>
`}(t);return s.html`
        <template
            role="row"
            class="${t=>"default"!==t.rowType?t.rowType:""}"
            :defaultCellItemTemplate="${i}"
            :defaultHeaderCellItemTemplate="${o}"
            ${(0,s.children)({property:"cellElements",filter:(0,s.elements)('[role="cell"],[role="gridcell"],[role="columnheader"],[role="rowheader"]')})}
        >
            <slot ${(0,s.slotted)("slottedCellElements")}></slot>
        </template>
    `},ii=(t,e)=>s.html`
        <template
            tabindex="-1"
            role="${t=>t.cellType&&"default"!==t.cellType?t.cellType:"gridcell"}"
            class="
            ${t=>"columnheader"===t.cellType?"column-header":"rowheader"===t.cellType?"row-header":""}
            "
        >
            <slot></slot>
        </template>
    `;function si(t){const e=t.parentElement;if(e)return e;{const e=t.getRootNode();if(e.host instanceof HTMLElement)return e.host}return null}function oi(t,e){let i=e;for(;null!==i;){if(i===t)return!0;i=si(i)}return!1}const ni=document.createElement("div");class ai{setProperty(t,e){s.DOM.queueUpdate((()=>this.target.setProperty(t,e)))}removeProperty(t){s.DOM.queueUpdate((()=>this.target.removeProperty(t)))}}class ri extends ai{constructor(){super();const t=new CSSStyleSheet;this.target=t.cssRules[t.insertRule(":root{}")].style,document.adoptedStyleSheets=[...document.adoptedStyleSheets,t]}}class li extends ai{constructor(){super(),this.style=document.createElement("style"),document.head.appendChild(this.style);const{sheet:t}=this.style;if(t){const e=t.insertRule(":root{}",t.cssRules.length);this.target=t.cssRules[e].style}}}class hi{constructor(t){this.store=new Map,this.target=null;const e=t.$fastController;this.style=document.createElement("style"),e.addStyles(this.style),s.Observable.getNotifier(e).subscribe(this,"isConnected"),this.handleChange(e,"isConnected")}targetChanged(){if(null!==this.target)for(const[t,e]of this.store.entries())this.target.setProperty(t,e)}setProperty(t,e){this.store.set(t,e),s.DOM.queueUpdate((()=>{null!==this.target&&this.target.setProperty(t,e)}))}removeProperty(t){this.store.delete(t),s.DOM.queueUpdate((()=>{null!==this.target&&this.target.removeProperty(t)}))}handleChange(t,e){const{sheet:i}=this.style;if(i){const t=i.insertRule(":host{}",i.cssRules.length);this.target=i.cssRules[t].style}else this.target=null}}(0,d.gn)([s.observable],hi.prototype,"target",void 0);class di{constructor(t){this.target=t.style}setProperty(t,e){s.DOM.queueUpdate((()=>this.target.setProperty(t,e)))}removeProperty(t){s.DOM.queueUpdate((()=>this.target.removeProperty(t)))}}class ci{setProperty(t,e){ci.properties[t]=e;for(const i of ci.roots.values())mi.getOrCreate(ci.normalizeRoot(i)).setProperty(t,e)}removeProperty(t){delete ci.properties[t];for(const e of ci.roots.values())mi.getOrCreate(ci.normalizeRoot(e)).removeProperty(t)}static registerRoot(t){const{roots:e}=ci;if(!e.has(t)){e.add(t);const i=mi.getOrCreate(this.normalizeRoot(t));for(const t in ci.properties)i.setProperty(t,ci.properties[t])}}static unregisterRoot(t){const{roots:e}=ci;if(e.has(t)){e.delete(t);const i=mi.getOrCreate(ci.normalizeRoot(t));for(const t in ci.properties)i.removeProperty(t)}}static normalizeRoot(t){return t===ni?document:t}}ci.roots=new Set,ci.properties={};const ui=new WeakMap,pi=s.DOM.supportsAdoptedStyleSheets?class extends ai{constructor(t){super();const e=new CSSStyleSheet;this.target=e.cssRules[e.insertRule(":host{}")].style,t.$fastController.addStyles(s.ElementStyles.create([e]))}}:hi,mi=Object.freeze({getOrCreate(t){if(ui.has(t))return ui.get(t);let e;return e=t===ni?new ci:t instanceof Document?s.DOM.supportsAdoptedStyleSheets?new ri:new li:t instanceof s.FASTElement?new pi(t):new di(t),ui.set(t,e),e}});class vi extends s.CSSDirective{constructor(t){super(),this.subscribers=new WeakMap,this._appliedTo=new Set,this.name=t.name,null!==t.cssCustomPropertyName&&(this.cssCustomProperty=`--${t.cssCustomPropertyName}`,this.cssVar=`var(${this.cssCustomProperty})`),this.id=vi.uniqueId(),vi.tokensById.set(this.id,this)}get appliedTo(){return[...this._appliedTo]}static from(t){return new vi({name:"string"==typeof t?t:t.name,cssCustomPropertyName:"string"==typeof t?t:void 0===t.cssCustomPropertyName?t.name:t.cssCustomPropertyName})}static isCSSDesignToken(t){return"string"==typeof t.cssCustomProperty}static isDerivedDesignTokenValue(t){return"function"==typeof t}static getTokenById(t){return vi.tokensById.get(t)}getOrCreateSubscriberSet(t=this){return this.subscribers.get(t)||this.subscribers.set(t,new Set)&&this.subscribers.get(t)}createCSS(){return this.cssVar||""}getValueFor(t){const e=Ci.getOrCreate(t).get(this);if(void 0!==e)return e;throw new Error(`Value could not be retrieved for token named "${this.name}". Ensure the value is set for ${t} or an ancestor of ${t}.`)}setValueFor(t,e){return this._appliedTo.add(t),e instanceof vi&&(e=this.alias(e)),Ci.getOrCreate(t).set(this,e),this}deleteValueFor(t){return this._appliedTo.delete(t),Ci.existsFor(t)&&Ci.getOrCreate(t).delete(this),this}withDefault(t){return this.setValueFor(ni,t),this}subscribe(t,e){const i=this.getOrCreateSubscriberSet(e);e&&!Ci.existsFor(e)&&Ci.getOrCreate(e),i.has(t)||i.add(t)}unsubscribe(t,e){const i=this.subscribers.get(e||this);i&&i.has(t)&&i.delete(t)}notify(t){const e=Object.freeze({token:this,target:t});this.subscribers.has(this)&&this.subscribers.get(this).forEach((t=>t.handleChange(e))),this.subscribers.has(t)&&this.subscribers.get(t).forEach((t=>t.handleChange(e)))}alias(t){return e=>t.getValueFor(e)}}vi.uniqueId=(()=>{let t=0;return()=>(t++,t.toString(16))})(),vi.tokensById=new Map;class gi{constructor(t,e,i){this.source=t,this.token=e,this.node=i,this.dependencies=new Set,this.observer=s.Observable.binding(t,this,!1),this.observer.handleChange=this.observer.call,this.handleChange()}disconnect(){this.observer.disconnect()}handleChange(){this.node.store.set(this.token,this.observer.observe(this.node.target,s.defaultExecutionContext))}}class bi{constructor(){this.values=new Map}set(t,e){this.values.get(t)!==e&&(this.values.set(t,e),s.Observable.getNotifier(this).notify(t.id))}get(t){return s.Observable.track(this,t.id),this.values.get(t)}delete(t){this.values.delete(t)}all(){return this.values.entries()}}const fi=new WeakMap,yi=new WeakMap;class Ci{constructor(t){this.target=t,this.store=new bi,this.children=[],this.assignedValues=new Map,this.reflecting=new Set,this.bindingObservers=new Map,this.tokenValueChangeHandler={handleChange:(t,e)=>{const i=vi.getTokenById(e);if(i&&(i.notify(this.target),vi.isCSSDesignToken(i))){const e=this.parent,s=this.isReflecting(i);if(e){const o=e.get(i),n=t.get(i);o===n||s?o===n&&s&&this.stopReflectToCSS(i):this.reflectToCSS(i)}else s||this.reflectToCSS(i)}}},fi.set(t,this),s.Observable.getNotifier(this.store).subscribe(this.tokenValueChangeHandler),t instanceof s.FASTElement?t.$fastController.addBehaviors([this]):t.isConnected&&this.bind()}static getOrCreate(t){return fi.get(t)||new Ci(t)}static existsFor(t){return fi.has(t)}static findParent(t){if(ni!==t.target){let e=si(t.target);for(;null!==e;){if(fi.has(e))return fi.get(e);e=si(e)}return Ci.getOrCreate(ni)}return null}static findClosestAssignedNode(t,e){let i=e;do{if(i.has(t))return i;i=i.parent?i.parent:i.target!==ni?Ci.getOrCreate(ni):null}while(null!==i);return null}get parent(){return yi.get(this)||null}has(t){return this.assignedValues.has(t)}get(t){const e=this.store.get(t);if(void 0!==e)return e;const i=this.getRaw(t);return void 0!==i?(this.hydrate(t,i),this.get(t)):void 0}getRaw(t){var e;return this.assignedValues.has(t)?this.assignedValues.get(t):null===(e=Ci.findClosestAssignedNode(t,this))||void 0===e?void 0:e.getRaw(t)}set(t,e){vi.isDerivedDesignTokenValue(this.assignedValues.get(t))&&this.tearDownBindingObserver(t),this.assignedValues.set(t,e),vi.isDerivedDesignTokenValue(e)?this.setupBindingObserver(t,e):this.store.set(t,e)}delete(t){this.assignedValues.delete(t),this.tearDownBindingObserver(t);const e=this.getRaw(t);e?this.hydrate(t,e):this.store.delete(t)}bind(){const t=Ci.findParent(this);t&&t.appendChild(this);for(const t of this.assignedValues.keys())t.notify(this.target)}unbind(){this.parent&&yi.get(this).removeChild(this)}appendChild(t){t.parent&&yi.get(t).removeChild(t);const e=this.children.filter((e=>t.contains(e)));yi.set(t,this),this.children.push(t),e.forEach((e=>t.appendChild(e))),s.Observable.getNotifier(this.store).subscribe(t);for(const[e,i]of this.store.all())t.hydrate(e,this.bindingObservers.has(e)?this.getRaw(e):i)}removeChild(t){const e=this.children.indexOf(t);return-1!==e&&this.children.splice(e,1),s.Observable.getNotifier(this.store).unsubscribe(t),t.parent===this&&yi.delete(t)}contains(t){return oi(this.target,t.target)}reflectToCSS(t){this.isReflecting(t)||(this.reflecting.add(t),Ci.cssCustomPropertyReflector.startReflection(t,this.target))}stopReflectToCSS(t){this.isReflecting(t)&&(this.reflecting.delete(t),Ci.cssCustomPropertyReflector.stopReflection(t,this.target))}isReflecting(t){return this.reflecting.has(t)}handleChange(t,e){const i=vi.getTokenById(e);i&&this.hydrate(i,this.getRaw(i))}hydrate(t,e){if(!this.has(t)){const i=this.bindingObservers.get(t);vi.isDerivedDesignTokenValue(e)?i?i.source!==e&&(this.tearDownBindingObserver(t),this.setupBindingObserver(t,e)):this.setupBindingObserver(t,e):(i&&this.tearDownBindingObserver(t),this.store.set(t,e))}}setupBindingObserver(t,e){const i=new gi(e,t,this);return this.bindingObservers.set(t,i),i}tearDownBindingObserver(t){return!!this.bindingObservers.has(t)&&(this.bindingObservers.get(t).disconnect(),this.bindingObservers.delete(t),!0)}}Ci.cssCustomPropertyReflector=new class{startReflection(t,e){t.subscribe(this,e),this.handleChange({token:t,target:e})}stopReflection(t,e){t.unsubscribe(this,e),this.remove(t,e)}handleChange(t){const{token:e,target:i}=t;this.add(e,i)}add(t,e){mi.getOrCreate(e).setProperty(t.cssCustomProperty,this.resolveCSSValue(Ci.getOrCreate(e).get(t)))}remove(t,e){mi.getOrCreate(e).removeProperty(t.cssCustomProperty)}resolveCSSValue(t){return t&&"function"==typeof t.createCSS?t.createCSS():t}},(0,d.gn)([s.observable],Ci.prototype,"children",void 0);const xi=Object.freeze({create:function(t){return vi.from(t)},notifyConnection:t=>!(!t.isConnected||!Ci.existsFor(t)||(Ci.getOrCreate(t).bind(),0)),notifyDisconnection:t=>!(t.isConnected||!Ci.existsFor(t)||(Ci.getOrCreate(t).unbind(),0)),registerRoot(t=ni){ci.registerRoot(t)},unregisterRoot(t=ni){ci.unregisterRoot(t)}}),wi=Object.freeze({definitionCallbackOnly:null,ignoreDuplicate:Symbol()}),$i=new Map,Ii=new Map;let ki=null;const Ei=y.createInterface((t=>t.cachedCallback((t=>(null===ki&&(ki=new Oi(null,t)),ki))))),Ti=Object.freeze({tagFor:t=>Ii.get(t),responsibleFor(t){const e=t.$$designSystem$$;return e||y.findResponsibleContainer(t).get(Ei)},getOrCreate(t){if(!t)return null===ki&&(ki=y.getOrCreateDOMContainer().get(Ei)),ki;const e=t.$$designSystem$$;if(e)return e;const i=y.getOrCreateDOMContainer(t);if(i.has(Ei,!1))return i.get(Ei);{const e=new Oi(t,i);return i.register(Y.instance(Ei,e)),e}}});class Oi{constructor(t,e){this.owner=t,this.container=e,this.designTokensInitialized=!1,this.prefix="fast",this.shadowRootMode=void 0,this.disambiguate=()=>wi.definitionCallbackOnly,null!==t&&(t.$$designSystem$$=this)}withPrefix(t){return this.prefix=t,this}withShadowRootMode(t){return this.shadowRootMode=t,this}withElementDisambiguation(t){return this.disambiguate=t,this}withDesignTokenRoot(t){return this.designTokenRoot=t,this}register(...t){const e=this.container,i=[],s=this.disambiguate,o=this.shadowRootMode,n={elementPrefix:this.prefix,tryDefineElement(t,n,a){const r=function(t,e,i){return"string"==typeof t?{name:t,type:e,callback:i}:t}(t,n,a),{name:l,callback:h,baseClass:d}=r;let{type:c}=r,u=l,p=$i.get(u),m=!0;for(;p;){const t=s(u,c,p);switch(t){case wi.ignoreDuplicate:return;case wi.definitionCallbackOnly:m=!1,p=void 0;break;default:u=t,p=$i.get(u)}}m&&((Ii.has(c)||c===rt)&&(c=class extends c{}),$i.set(u,c),Ii.set(c,u),d&&Ii.set(d,u)),i.push(new Di(e,u,c,o,h,m))}};this.designTokensInitialized||(this.designTokensInitialized=!0,null!==this.designTokenRoot&&xi.registerRoot(this.designTokenRoot)),e.registerWithContext(n,...t);for(const t of i)t.callback(t),t.willDefine&&null!==t.definition&&t.definition.define();return this}}class Di{constructor(t,e,i,s,o,n){this.container=t,this.name=e,this.type=i,this.shadowRootMode=s,this.callback=o,this.willDefine=n,this.definition=null}definePresentation(t){nt.define(this.name,t,this.container)}defineElement(t){this.definition=new s.FASTElementDefinition(this.type,Object.assign(Object.assign({},t),{name:this.name}))}tagFor(t){return Ti.tagFor(t)}}const Ri=(t,e)=>s.html`
    <div class="positioning-region" part="positioning-region">
        ${(0,s.when)((t=>t.modal),s.html`
                <div
                    class="overlay"
                    part="overlay"
                    role="presentation"
                    @click="${t=>t.dismiss()}"
                ></div>
            `)}
        <div
            role="dialog"
            tabindex="-1"
            class="control"
            part="control"
            aria-modal="${t=>t.modal}"
            aria-describedby="${t=>t.ariaDescribedby}"
            aria-labelledby="${t=>t.ariaLabelledby}"
            aria-label="${t=>t.ariaLabel}"
            ${(0,s.ref)("dialog")}
        >
            <slot></slot>
        </div>
    </div>
`;var Si=["input","select","textarea","a[href]","button","[tabindex]:not(slot)","audio[controls]","video[controls]",'[contenteditable]:not([contenteditable="false"])',"details>summary:first-of-type","details"],Ai=Si.join(","),Fi="undefined"==typeof Element,Li=Fi?function(){}:Element.prototype.matches||Element.prototype.msMatchesSelector||Element.prototype.webkitMatchesSelector,Pi=!Fi&&Element.prototype.getRootNode?function(t){return t.getRootNode()}:function(t){return t.ownerDocument},Mi=function(t){return"INPUT"===t.tagName},Hi=function(t){var e=t.getBoundingClientRect(),i=e.width,s=e.height;return 0===i&&0===s},zi=function(t,e){return!(e.disabled||function(t){return Mi(t)&&"hidden"===t.type}(e)||function(t,e){var i=e.displayCheck,s=e.getShadowRoot;if("hidden"===getComputedStyle(t).visibility)return!0;var o=Li.call(t,"details>summary:first-of-type")?t.parentElement:t;if(Li.call(o,"details:not([open]) *"))return!0;var n=Pi(t).host,a=(null==n?void 0:n.ownerDocument.contains(n))||t.ownerDocument.contains(t);if(i&&"full"!==i){if("non-zero-area"===i)return Hi(t)}else{if("function"==typeof s){for(var r=t;t;){var l=t.parentElement,h=Pi(t);if(l&&!l.shadowRoot&&!0===s(l))return Hi(t);t=t.assignedSlot?t.assignedSlot:l||h===t.ownerDocument?l:h.host}t=r}if(a)return!t.getClientRects().length}return!1}(e,t)||function(t){return"DETAILS"===t.tagName&&Array.prototype.slice.apply(t.children).some((function(t){return"SUMMARY"===t.tagName}))}(e)||function(t){if(/^(INPUT|BUTTON|SELECT|TEXTAREA)$/.test(t.tagName))for(var e=t.parentElement;e;){if("FIELDSET"===e.tagName&&e.disabled){for(var i=0;i<e.children.length;i++){var s=e.children.item(i);if("LEGEND"===s.tagName)return!!Li.call(e,"fieldset[disabled] *")||!s.contains(t)}return!0}e=e.parentElement}return!1}(e))},Vi=function(t,e){return!(function(t){return function(t){return Mi(t)&&"radio"===t.type}(t)&&!function(t){if(!t.name)return!0;var e,i=t.form||Pi(t),s=function(t){return i.querySelectorAll('input[type="radio"][name="'+t+'"]')};if("undefined"!=typeof window&&void 0!==window.CSS&&"function"==typeof window.CSS.escape)e=s(window.CSS.escape(t.name));else try{e=s(t.name)}catch(t){return console.error("Looks like you have a radio button with a name attribute containing invalid CSS selector characters and need the CSS.escape polyfill: %s",t.message),!1}var o=function(t,e){for(var i=0;i<t.length;i++)if(t[i].checked&&t[i].form===e)return t[i]}(e,t.form);return!o||o===t}(t)}(e)||function(t,e){return t.tabIndex<0&&(e||/^(AUDIO|VIDEO|DETAILS)$/.test(t.tagName)||t.isContentEditable)&&isNaN(parseInt(t.getAttribute("tabindex"),10))?0:t.tabIndex}(e)<0||!zi(t,e))},Ni=function(t,e){if(e=e||{},!t)throw new Error("No node provided");return!1!==Li.call(t,Ai)&&Vi(e,t)},Bi=Si.concat("iframe").join(","),qi=function(t,e){if(e=e||{},!t)throw new Error("No node provided");return!1!==Li.call(t,Bi)&&zi(e,t)};class Ui extends rt{constructor(){super(...arguments),this.modal=!0,this.hidden=!1,this.trapFocus=!0,this.trapFocusChanged=()=>{this.$fastController.isConnected&&this.updateTrapFocus()},this.isTrappingFocus=!1,this.handleDocumentKeydown=t=>{if(!t.defaultPrevented&&!this.hidden)switch(t.key){case yt:this.dismiss(),t.preventDefault();break;case $t:this.handleTabKeyDown(t)}},this.handleDocumentFocus=t=>{!t.defaultPrevented&&this.shouldForceFocus(t.target)&&(this.focusFirstElement(),t.preventDefault())},this.handleTabKeyDown=t=>{if(!this.trapFocus||this.hidden)return;const e=this.getTabQueueBounds();return 0!==e.length?1===e.length?(e[0].focus(),void t.preventDefault()):void(t.shiftKey&&t.target===e[0]?(e[e.length-1].focus(),t.preventDefault()):t.shiftKey||t.target!==e[e.length-1]||(e[0].focus(),t.preventDefault())):void 0},this.getTabQueueBounds=()=>Ui.reduceTabbableItems([],this),this.focusFirstElement=()=>{const t=this.getTabQueueBounds();t.length>0?t[0].focus():this.dialog instanceof HTMLElement&&this.dialog.focus()},this.shouldForceFocus=t=>this.isTrappingFocus&&!this.contains(t),this.shouldTrapFocus=()=>this.trapFocus&&!this.hidden,this.updateTrapFocus=t=>{const e=void 0===t?this.shouldTrapFocus():t;e&&!this.isTrappingFocus?(this.isTrappingFocus=!0,document.addEventListener("focusin",this.handleDocumentFocus),s.DOM.queueUpdate((()=>{this.shouldForceFocus(document.activeElement)&&this.focusFirstElement()}))):!e&&this.isTrappingFocus&&(this.isTrappingFocus=!1,document.removeEventListener("focusin",this.handleDocumentFocus))}}dismiss(){this.$emit("dismiss"),this.$emit("cancel")}show(){this.hidden=!1}hide(){this.hidden=!0,this.$emit("close")}connectedCallback(){super.connectedCallback(),document.addEventListener("keydown",this.handleDocumentKeydown),this.notifier=s.Observable.getNotifier(this),this.notifier.subscribe(this,"hidden"),this.updateTrapFocus()}disconnectedCallback(){super.disconnectedCallback(),document.removeEventListener("keydown",this.handleDocumentKeydown),this.updateTrapFocus(!1),this.notifier.unsubscribe(this,"hidden")}handleChange(t,e){"hidden"===e&&this.updateTrapFocus()}static reduceTabbableItems(t,e){return"-1"===e.getAttribute("tabindex")?t:Ni(e)||Ui.isFocusableFastElement(e)&&Ui.hasTabbableShadow(e)?(t.push(e),t):e.childElementCount?t.concat(Array.from(e.children).reduce(Ui.reduceTabbableItems,[])):t}static isFocusableFastElement(t){var e,i;return!!(null===(i=null===(e=t.$fastController)||void 0===e?void 0:e.definition.shadowOptions)||void 0===i?void 0:i.delegatesFocus)}static hasTabbableShadow(t){var e,i;return Array.from(null!==(i=null===(e=t.shadowRoot)||void 0===e?void 0:e.querySelectorAll("*"))&&void 0!==i?i:[]).some((t=>Ni(t)))}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],Ui.prototype,"modal",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Ui.prototype,"hidden",void 0),(0,d.gn)([(0,s.attr)({attribute:"trap-focus",mode:"boolean"})],Ui.prototype,"trapFocus",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-describedby"})],Ui.prototype,"ariaDescribedby",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-labelledby"})],Ui.prototype,"ariaLabelledby",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-label"})],Ui.prototype,"ariaLabel",void 0);const ji=new MutationObserver((t=>{for(const e of t)_i.getOrCreateFor(e.target).notify(e.attributeName)}));class _i extends s.SubscriberSet{constructor(t){super(t),this.watchedAttributes=new Set,_i.subscriberCache.set(t,this)}subscribe(t){super.subscribe(t),this.watchedAttributes.has(t.attributes)||(this.watchedAttributes.add(t.attributes),this.observe())}unsubscribe(t){super.unsubscribe(t),this.watchedAttributes.has(t.attributes)&&(this.watchedAttributes.delete(t.attributes),this.observe())}static getOrCreateFor(t){return this.subscriberCache.get(t)||new _i(t)}observe(){const t=[];for(const e of this.watchedAttributes.values())for(let i=0;i<e.length;i++)t.push(e[i]);ji.observe(this.source,{attributeFilter:t})}}_i.subscriberCache=new WeakMap;class Ki{constructor(t,e){this.target=t,this.attributes=Object.freeze(e)}bind(t){if(_i.getOrCreateFor(t).subscribe(this),t.hasAttributes())for(let e=0;e<t.attributes.length;e++)this.handleChange(t,t.attributes[e].name)}unbind(t){_i.getOrCreateFor(t).unsubscribe(this)}handleChange(t,e){this.attributes.includes(e)&&s.DOM.setAttribute(this.target,e,t.getAttribute(e))}}function Wi(...t){return new s.AttachedBehaviorHTMLDirective("fast-reflect-attr",Ki,t)}const Gi=(t,e)=>s.html`
    <details class="disclosure" ${(0,s.ref)("details")}>
        <summary
            class="invoker"
            role="button"
            aria-controls="disclosure-content"
            aria-expanded="${t=>t.expanded}"
        >
            <slot name="start"></slot>
            <slot name="title">${t=>t.title}</slot>
            <slot name="end"></slot>
        </summary>
        <div id="disclosure-content"><slot></slot></div>
    </details>
`;class Yi extends rt{connectedCallback(){super.connectedCallback(),this.setup()}disconnectedCallback(){super.disconnectedCallback(),this.details.removeEventListener("toggle",this.onToggle)}show(){this.details.open=!0}hide(){this.details.open=!1}toggle(){this.details.open=!this.details.open}setup(){this.onToggle=this.onToggle.bind(this),this.details.addEventListener("toggle",this.onToggle),this.expanded&&this.show()}onToggle(){this.expanded=this.details.open,this.$emit("toggle")}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],Yi.prototype,"expanded",void 0),(0,d.gn)([s.attr],Yi.prototype,"title",void 0);const Xi=(t,e)=>s.html`
    <template role="${t=>t.role}" aria-orientation="${t=>t.orientation}"></template>
`;var Qi=i(2899);const Zi={separator:"separator",presentation:"presentation"};class Ji extends rt{constructor(){super(...arguments),this.role=Zi.separator,this.orientation=Qi.i.horizontal}}(0,d.gn)([s.attr],Ji.prototype,"role",void 0),(0,d.gn)([s.attr],Ji.prototype,"orientation",void 0);const ts={next:"next",previous:"previous"},es=(t,e)=>s.html`
    <template
        role="button"
        aria-disabled="${t=>!!t.disabled||void 0}"
        tabindex="${t=>t.hiddenFromAT?-1:0}"
        class="${t=>t.direction} ${t=>t.disabled?"disabled":""}"
        @keyup="${(t,e)=>t.keyupHandler(e.event)}"
    >
        ${(0,s.when)((t=>t.direction===ts.next),s.html`
                <span part="next" class="next">
                    <slot name="next">
                        ${e.next||""}
                    </slot>
                </span>
            `)}
        ${(0,s.when)((t=>t.direction===ts.previous),s.html`
                <span part="previous" class="previous">
                    <slot name="previous">
                        ${e.previous||""}
                    </slot>
                </span>
            `)}
    </template>
`;class is extends rt{constructor(){super(...arguments),this.hiddenFromAT=!0,this.direction=ts.next}keyupHandler(t){if(!this.hiddenFromAT){const e=t.key;"Enter"!==e&&"Space"!==e||this.$emit("click",t),"Escape"===e&&this.blur()}}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],is.prototype,"disabled",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-hidden",converter:s.booleanConverter})],is.prototype,"hiddenFromAT",void 0),(0,d.gn)([s.attr],is.prototype,"direction",void 0);const ss=(t,e)=>s.html`
    <template
        aria-checked="${t=>t.ariaChecked}"
        aria-disabled="${t=>t.ariaDisabled}"
        aria-posinset="${t=>t.ariaPosInSet}"
        aria-selected="${t=>t.ariaSelected}"
        aria-setsize="${t=>t.ariaSetSize}"
        class="${t=>[t.checked&&"checked",t.selected&&"selected",t.disabled&&"disabled"].filter(Boolean).join(" ")}"
        role="option"
    >
        ${a(t,e)}
        <span class="content" part="content">
            <slot ${(0,s.slotted)("content")}></slot>
        </span>
        ${n(t,e)}
    </template>
`;class os extends _e{constructor(){super(...arguments),this.activeIndex=-1,this.rangeStartIndex=-1}get activeOption(){return this.options[this.activeIndex]}get checkedOptions(){var t;return null===(t=this.options)||void 0===t?void 0:t.filter((t=>t.checked))}get firstSelectedOptionIndex(){return this.options.indexOf(this.firstSelectedOption)}activeIndexChanged(t,e){var i,s;this.ariaActiveDescendant=null!==(s=null===(i=this.options[e])||void 0===i?void 0:i.id)&&void 0!==s?s:"",this.focusAndScrollOptionIntoView()}checkActiveIndex(){if(!this.multiple)return;const t=this.activeOption;t&&(t.checked=!0)}checkFirstOption(t=!1){t?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex+1),this.options.forEach(((t,e)=>{t.checked=Et(e,this.rangeStartIndex)}))):this.uncheckAllOptions(),this.activeIndex=0,this.checkActiveIndex()}checkLastOption(t=!1){t?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),this.options.forEach(((t,e)=>{t.checked=Et(e,this.rangeStartIndex,this.options.length)}))):this.uncheckAllOptions(),this.activeIndex=this.options.length-1,this.checkActiveIndex()}connectedCallback(){super.connectedCallback(),this.addEventListener("focusout",this.focusoutHandler)}disconnectedCallback(){this.removeEventListener("focusout",this.focusoutHandler),super.disconnectedCallback()}checkNextOption(t=!1){t?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),this.options.forEach(((t,e)=>{t.checked=Et(e,this.rangeStartIndex,this.activeIndex+1)}))):this.uncheckAllOptions(),this.activeIndex+=this.activeIndex<this.options.length-1?1:0,this.checkActiveIndex()}checkPreviousOption(t=!1){t?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),1===this.checkedOptions.length&&(this.rangeStartIndex+=1),this.options.forEach(((t,e)=>{t.checked=Et(e,this.activeIndex,this.rangeStartIndex)}))):this.uncheckAllOptions(),this.activeIndex-=this.activeIndex>0?1:0,this.checkActiveIndex()}clickHandler(t){var e;if(!this.multiple)return super.clickHandler(t);const i=null===(e=t.target)||void 0===e?void 0:e.closest("[role=option]");return i&&!i.disabled?(this.uncheckAllOptions(),this.activeIndex=this.options.indexOf(i),this.checkActiveIndex(),this.toggleSelectedForAllCheckedOptions(),!0):void 0}focusAndScrollOptionIntoView(){super.focusAndScrollOptionIntoView(this.activeOption)}focusinHandler(t){if(!this.multiple)return super.focusinHandler(t);this.shouldSkipFocus||t.target!==t.currentTarget||(this.uncheckAllOptions(),-1===this.activeIndex&&(this.activeIndex=-1!==this.firstSelectedOptionIndex?this.firstSelectedOptionIndex:0),this.checkActiveIndex(),this.setSelectedOptions(),this.focusAndScrollOptionIntoView()),this.shouldSkipFocus=!1}focusoutHandler(t){this.multiple&&this.uncheckAllOptions()}keydownHandler(t){if(!this.multiple)return super.keydownHandler(t);if(this.disabled)return!0;const{key:e,shiftKey:i}=t;switch(this.shouldSkipFocus=!1,e){case Ct:return void this.checkFirstOption(i);case mt:return void this.checkNextOption(i);case bt:return void this.checkPreviousOption(i);case xt:return void this.checkLastOption(i);case $t:return this.focusAndScrollOptionIntoView(),!0;case yt:return this.uncheckAllOptions(),this.checkActiveIndex(),!0;case wt:if(t.preventDefault(),this.typeAheadExpired)return void this.toggleSelectedForAllCheckedOptions();default:return 1===e.length&&this.handleTypeAhead(`${e}`),!0}}mousedownHandler(t){if(t.offsetX>=0&&t.offsetX<=this.scrollWidth)return super.mousedownHandler(t)}multipleChanged(t,e){var i;this.ariaMultiSelectable=e?"true":null,null===(i=this.options)||void 0===i||i.forEach((t=>{t.checked=!e&&void 0})),this.setSelectedOptions()}setSelectedOptions(){this.multiple?this.$fastController.isConnected&&this.options&&(this.selectedOptions=this.options.filter((t=>t.selected)),this.focusAndScrollOptionIntoView()):super.setSelectedOptions()}sizeChanged(t,e){var i;const o=Math.max(0,parseInt(null!==(i=null==e?void 0:e.toFixed())&&void 0!==i?i:"",10));o!==e&&s.DOM.queueUpdate((()=>{this.size=o}))}toggleSelectedForAllCheckedOptions(){const t=this.checkedOptions.filter((t=>!t.disabled)),e=!t.every((t=>t.selected));t.forEach((t=>t.selected=e)),this.selectedIndex=this.options.indexOf(t[t.length-1]),this.setSelectedOptions()}typeaheadBufferChanged(t,e){if(this.multiple){if(this.$fastController.isConnected){const t=this.getTypeaheadMatches(),e=this.options.indexOf(t[0]);e>-1&&(this.activeIndex=e,this.uncheckAllOptions(),this.checkActiveIndex()),this.typeAheadExpired=!1}}else super.typeaheadBufferChanged(t,e)}uncheckAllOptions(t=!1){this.options.forEach((t=>t.checked=!this.multiple&&void 0)),t||(this.rangeStartIndex=-1)}}(0,d.gn)([s.observable],os.prototype,"activeIndex",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],os.prototype,"multiple",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],os.prototype,"size",void 0);const ns=(t,e)=>s.html`
    <template
        aria-activedescendant="${t=>t.ariaActiveDescendant}"
        aria-multiselectable="${t=>t.ariaMultiSelectable}"
        class="listbox"
        role="listbox"
        tabindex="${t=>t.disabled?null:"0"}"
        @click="${(t,e)=>t.clickHandler(e.event)}"
        @focusin="${(t,e)=>t.focusinHandler(e.event)}"
        @keydown="${(t,e)=>t.keydownHandler(e.event)}"
        @mousedown="${(t,e)=>t.mousedownHandler(e.event)}"
    >
        <slot
            ${(0,s.slotted)({filter:os.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
        ></slot>
    </template>
`;class as extends rt{constructor(){super(...arguments),this.optionElements=[]}menuElementsChanged(){this.updateOptions()}headerElementsChanged(){this.updateOptions()}footerElementsChanged(){this.updateOptions()}updateOptions(){this.optionElements.splice(0,this.optionElements.length),this.addSlottedListItems(this.headerElements),this.addSlottedListItems(this.menuElements),this.addSlottedListItems(this.footerElements),this.$emit("optionsupdated",{bubbles:!1})}addSlottedListItems(t){void 0!==t&&t.forEach((t=>{1===t.nodeType&&"listitem"===t.getAttribute("role")&&(t.id=t.id||Ne("option-"),this.optionElements.push(t))}))}}(0,d.gn)([s.observable],as.prototype,"menuElements",void 0),(0,d.gn)([s.observable],as.prototype,"headerElements",void 0),(0,d.gn)([s.observable],as.prototype,"footerElements",void 0),(0,d.gn)([s.observable],as.prototype,"suggestionsAvailableText",void 0);const rs=s.html`
    <template>
        ${t=>t.value}
    </template>
`;class ls extends rt{contentsTemplateChanged(){this.$fastController.isConnected&&this.updateView()}connectedCallback(){super.connectedCallback(),this.updateView()}disconnectedCallback(){super.disconnectedCallback(),this.disconnectView()}handleClick(t){return t.defaultPrevented||this.handleInvoked(),!1}handleInvoked(){this.$emit("pickeroptioninvoked")}updateView(){var t,e;this.disconnectView(),this.customView=null!==(e=null===(t=this.contentsTemplate)||void 0===t?void 0:t.render(this,this))&&void 0!==e?e:rs.render(this,this)}disconnectView(){var t;null===(t=this.customView)||void 0===t||t.dispose(),this.customView=void 0}}(0,d.gn)([(0,s.attr)({attribute:"value"})],ls.prototype,"value",void 0),(0,d.gn)([s.observable],ls.prototype,"contentsTemplate",void 0);class hs extends rt{}const ds=s.html`
    <template>
        ${t=>t.value}
    </template>
`;class cs extends rt{contentsTemplateChanged(){this.$fastController.isConnected&&this.updateView()}connectedCallback(){super.connectedCallback(),this.updateView()}disconnectedCallback(){this.disconnectView(),super.disconnectedCallback()}handleKeyDown(t){return!(t.defaultPrevented||t.key===ft&&(this.handleInvoke(),1))}handleClick(t){return t.defaultPrevented||this.handleInvoke(),!1}handleInvoke(){this.$emit("pickeriteminvoked")}updateView(){var t,e;this.disconnectView(),this.customView=null!==(e=null===(t=this.contentsTemplate)||void 0===t?void 0:t.render(this,this))&&void 0!==e?e:ds.render(this,this)}disconnectView(){var t;null===(t=this.customView)||void 0===t||t.dispose(),this.customView=void 0}}(0,d.gn)([(0,s.attr)({attribute:"value"})],cs.prototype,"value",void 0),(0,d.gn)([s.observable],cs.prototype,"contentsTemplate",void 0);const us=(t,e)=>{const i=t.tagFor(qt),o=t.tagFor(as),n=t.tagFor(hs),a=t.tagFor(hs),r=function(t){const e=t.tagFor(cs);return s.html`
    <${e}
        value="${t=>t}"
        :contentsTemplate="${(t,e)=>e.parent.listItemContentsTemplate}"
    >
    </${e}>
    `}(t),l=function(t){const e=t.tagFor(ls);return s.html`
    <${e}
        value="${t=>t}"
        :contentsTemplate="${(t,e)=>e.parent.menuOptionContentsTemplate}"
    >
    </${e}>
    `}(t);return s.html`
        <template
            :selectedListTag="${()=>n}"
            :menuTag="${()=>o}"
            :defaultListItemTemplate="${r}"
            :defaultMenuOptionTemplate="${l}"
            @focusin="${(t,e)=>t.handleFocusIn(e.event)}"
            @focusout="${(t,e)=>t.handleFocusOut(e.event)}"
            @keydown="${(t,e)=>t.handleKeyDown(e.event)}"
            @pickeriteminvoked="${(t,e)=>t.handleItemInvoke(e.event)}"
            @pickeroptioninvoked="${(t,e)=>t.handleOptionInvoke(e.event)}"
        >
            <slot name="list-region"></slot>

            ${(0,s.when)((t=>t.flyoutOpen),s.html`
                <${i}
                    class="region"
                    part="region"
                    auto-update-mode="${t=>t.menuConfig.autoUpdateMode}"
                    fixed-placement="${t=>t.menuConfig.fixedPlacement}"
                    vertical-positioning-mode="${t=>t.menuConfig.verticalPositioningMode}"
                    vertical-default-position="${t=>t.menuConfig.verticalDefaultPosition}"
                    vertical-scaling="${t=>t.menuConfig.verticalScaling}"
                    vertical-inset="${t=>t.menuConfig.verticalInset}"
                    vertical-viewport-lock="${t=>t.menuConfig.verticalViewportLock}"
                    horizontal-positioning-mode="${t=>t.menuConfig.horizontalPositioningMode}"
                    horizontal-default-position="${t=>t.menuConfig.horizontalDefaultPosition}"
                    horizontal-scaling="${t=>t.menuConfig.horizontalScaling}"
                    horizontal-inset="${t=>t.menuConfig.horizontalInset}"
                    horizontal-viewport-lock="${t=>t.menuConfig.horizontalViewportLock}"
                    @loaded="${(t,e)=>t.handleRegionLoaded(e.event)}"
                    ${(0,s.ref)("region")}
                >
                    ${(0,s.when)((t=>!t.showNoOptions&&!t.showLoading),s.html`
                            <slot name="menu-region"></slot>
                        `)}
                    ${(0,s.when)((t=>t.showNoOptions&&!t.showLoading),s.html`
                            <div class="no-options-display" part="no-options-display">
                                <slot name="no-options-region">
                                    ${t=>t.noSuggestionsText}
                                </slot>
                            </div>
                        `)}
                    ${(0,s.when)((t=>t.showLoading),s.html`
                            <div class="loading-display" part="loading-display">
                                <slot name="loading-region">
                                    <${a}
                                        part="loading-progress"
                                        class="loading-progress
                                        slot="loading-region"
                                    ></${a}>
                                        ${t=>t.loadingText}
                                </slot>
                            </div>
                        `)}
                </${i}>
            `)}
        </template>
    `};class ps extends rt{}class ms extends(he(ps)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}const vs=s.html`
    <input
        slot="input-region"
        role="combobox"
        type="text"
        autocapitalize="off"
        autocomplete="off"
        haspopup="list"
        aria-label="${t=>t.label}"
        aria-labelledby="${t=>t.labelledBy}"
        placeholder="${t=>t.placeholder}"
        ${(0,s.ref)("inputElement")}
    ></input>
`;class gs extends ms{constructor(){super(...arguments),this.selection="",this.filterSelected=!0,this.filterQuery=!0,this.noSuggestionsText="No suggestions available",this.suggestionsAvailableText="Suggestions available",this.loadingText="Loading suggestions",this.menuPlacement="bottom-fill",this.showLoading=!1,this.optionsList=[],this.filteredOptionsList=[],this.flyoutOpen=!1,this.menuFocusIndex=-1,this.showNoOptions=!1,this.selectedItems=[],this.inputElementView=null,this.handleTextInput=t=>{this.query=this.inputElement.value},this.handleInputClick=t=>{t.preventDefault(),this.toggleFlyout(!0)},this.setRegionProps=()=>{this.flyoutOpen&&(null!==this.region&&void 0!==this.region?this.region.anchorElement=this.inputElement:s.DOM.queueUpdate(this.setRegionProps))},this.configLookup={top:jt,bottom:_t,tallest:Kt,"top-fill":Wt,"bottom-fill":Gt,"tallest-fill":Yt}}selectionChanged(){this.$fastController.isConnected&&(this.handleSelectionChange(),this.proxy instanceof HTMLInputElement&&(this.proxy.value=this.selection,this.validate()))}optionsChanged(){this.optionsList=this.options.split(",").map((t=>t.trim())).filter((t=>""!==t))}menuPlacementChanged(){this.$fastController.isConnected&&this.updateMenuConfig()}showLoadingChanged(){this.$fastController.isConnected&&s.DOM.queueUpdate((()=>{this.setFocusedOption(0)}))}listItemTemplateChanged(){this.updateListItemTemplate()}defaultListItemTemplateChanged(){this.updateListItemTemplate()}menuOptionTemplateChanged(){this.updateOptionTemplate()}defaultMenuOptionTemplateChanged(){this.updateOptionTemplate()}optionsListChanged(){this.updateFilteredOptions()}queryChanged(){this.$fastController.isConnected&&(this.inputElement.value!==this.query&&(this.inputElement.value=this.query),this.updateFilteredOptions(),this.$emit("querychange",{bubbles:!1}))}filteredOptionsListChanged(){this.$fastController.isConnected&&(this.showNoOptions=0===this.filteredOptionsList.length&&0===this.menuElement.querySelectorAll('[role="listitem"]').length,this.setFocusedOption(this.showNoOptions?-1:0))}flyoutOpenChanged(){this.flyoutOpen?(s.DOM.queueUpdate(this.setRegionProps),this.$emit("menuopening",{bubbles:!1})):this.$emit("menuclosing",{bubbles:!1})}showNoOptionsChanged(){this.$fastController.isConnected&&s.DOM.queueUpdate((()=>{this.setFocusedOption(0)}))}connectedCallback(){super.connectedCallback(),this.listElement=document.createElement(this.selectedListTag),this.appendChild(this.listElement),this.itemsPlaceholderElement=document.createComment(""),this.listElement.append(this.itemsPlaceholderElement),this.inputElementView=vs.render(this,this.listElement);const t=this.menuTag.toUpperCase();this.menuElement=Array.from(this.children).find((e=>e.tagName===t)),void 0===this.menuElement&&(this.menuElement=document.createElement(this.menuTag),this.appendChild(this.menuElement)),""===this.menuElement.id&&(this.menuElement.id=Ne("listbox-")),this.menuId=this.menuElement.id,this.optionsPlaceholder=document.createComment(""),this.menuElement.append(this.optionsPlaceholder),this.updateMenuConfig(),s.DOM.queueUpdate((()=>this.initialize()))}disconnectedCallback(){super.disconnectedCallback(),this.toggleFlyout(!1),this.inputElement.removeEventListener("input",this.handleTextInput),this.inputElement.removeEventListener("click",this.handleInputClick),null!==this.inputElementView&&(this.inputElementView.dispose(),this.inputElementView=null)}focus(){this.inputElement.focus()}initialize(){this.updateListItemTemplate(),this.updateOptionTemplate(),this.itemsRepeatBehavior=new s.RepeatDirective((t=>t.selectedItems),(t=>t.activeListItemTemplate),{positioning:!0}).createBehavior(this.itemsPlaceholderElement),this.inputElement.addEventListener("input",this.handleTextInput),this.inputElement.addEventListener("click",this.handleInputClick),this.$fastController.addBehaviors([this.itemsRepeatBehavior]),this.menuElement.suggestionsAvailableText=this.suggestionsAvailableText,this.menuElement.addEventListener("optionsupdated",this.handleMenuOptionsUpdated),this.optionsRepeatBehavior=new s.RepeatDirective((t=>t.filteredOptionsList),(t=>t.activeMenuOptionTemplate),{positioning:!0}).createBehavior(this.optionsPlaceholder),this.$fastController.addBehaviors([this.optionsRepeatBehavior]),this.handleSelectionChange()}toggleFlyout(t){if(this.flyoutOpen!==t){if(t&&document.activeElement===this.inputElement)return this.flyoutOpen=t,void s.DOM.queueUpdate((()=>{void 0!==this.menuElement?this.setFocusedOption(0):this.disableMenu()}));this.flyoutOpen=!1,this.disableMenu()}}handleMenuOptionsUpdated(t){t.preventDefault(),this.flyoutOpen&&this.setFocusedOption(0)}handleKeyDown(t){if(t.defaultPrevented)return!1;switch(t.key){case mt:if(this.flyoutOpen){const t=this.flyoutOpen?Math.min(this.menuFocusIndex+1,this.menuElement.optionElements.length-1):0;this.setFocusedOption(t)}else this.toggleFlyout(!0);return!1;case bt:if(this.flyoutOpen){const t=this.flyoutOpen?Math.max(this.menuFocusIndex-1,0):0;this.setFocusedOption(t)}else this.toggleFlyout(!0);return!1;case yt:return this.toggleFlyout(!1),!1;case ft:return-1!==this.menuFocusIndex&&this.menuElement.optionElements.length>this.menuFocusIndex&&this.menuElement.optionElements[this.menuFocusIndex].click(),!1;case gt:return document.activeElement===this.inputElement||(this.incrementFocusedItem(1),!1);case vt:return 0!==this.inputElement.selectionStart||(this.incrementFocusedItem(-1),!1);case"Delete":case"Backspace":{if(null===document.activeElement)return!0;if(document.activeElement===this.inputElement)return 0!==this.inputElement.selectionStart||(this.selection=this.selectedItems.slice(0,this.selectedItems.length-1).toString(),this.toggleFlyout(!1),!1);const t=Array.from(this.listElement.children),e=t.indexOf(document.activeElement);return!(e>-1&&(this.selection=this.selectedItems.splice(e,1).toString(),s.DOM.queueUpdate((()=>{t[Math.min(t.length,e)].focus()})),1))}}return this.toggleFlyout(!0),!0}handleFocusIn(t){return!1}handleFocusOut(t){return void 0!==this.menuElement&&this.menuElement.contains(t.relatedTarget)||this.toggleFlyout(!1),!1}handleSelectionChange(){this.selectedItems.toString()!==this.selection&&(this.selectedItems=""===this.selection?[]:this.selection.split(","),this.updateFilteredOptions(),s.DOM.queueUpdate((()=>{this.checkMaxItems()})),this.$emit("selectionchange",{bubbles:!1}))}handleRegionLoaded(t){s.DOM.queueUpdate((()=>{this.setFocusedOption(0),this.$emit("menuloaded",{bubbles:!1})}))}checkMaxItems(){if(void 0!==this.inputElement)if(void 0!==this.maxSelected&&this.selectedItems.length>=this.maxSelected){if(document.activeElement===this.inputElement){const t=Array.from(this.listElement.querySelectorAll("[role='listitem']"));t[t.length-1].focus()}this.inputElement.hidden=!0}else this.inputElement.hidden=!1}handleItemInvoke(t){if(t.defaultPrevented)return!1;if(t.target instanceof cs){const e=Array.from(this.listElement.querySelectorAll("[role='listitem']")).indexOf(t.target);if(-1!==e){const t=this.selectedItems.slice();t.splice(e,1),this.selection=t.toString(),s.DOM.queueUpdate((()=>this.incrementFocusedItem(0)))}return!1}return!0}handleOptionInvoke(t){return!(t.defaultPrevented||t.target instanceof ls&&(void 0!==t.target.value&&(this.selection=`${this.selection}${""===this.selection?"":","}${t.target.value}`),this.inputElement.value="",this.query="",this.inputElement.focus(),this.toggleFlyout(!1),1))}incrementFocusedItem(t){if(0===this.selectedItems.length)return void this.inputElement.focus();const e=Array.from(this.listElement.querySelectorAll("[role='listitem']"));if(null!==document.activeElement){let i=e.indexOf(document.activeElement);-1===i&&(i=e.length);const s=Math.min(e.length,Math.max(0,i+t));s===e.length?void 0!==this.maxSelected&&this.selectedItems.length>=this.maxSelected?e[s-1].focus():this.inputElement.focus():e[s].focus()}}disableMenu(){var t,e,i;this.menuFocusIndex=-1,this.menuFocusOptionId=void 0,null===(t=this.inputElement)||void 0===t||t.removeAttribute("aria-activedescendant"),null===(e=this.inputElement)||void 0===e||e.removeAttribute("aria-owns"),null===(i=this.inputElement)||void 0===i||i.removeAttribute("aria-expanded")}setFocusedOption(t){if(!this.flyoutOpen||-1===t||this.showNoOptions||this.showLoading)return void this.disableMenu();if(0===this.menuElement.optionElements.length)return;this.menuElement.optionElements.forEach((t=>{t.setAttribute("aria-selected","false")})),this.menuFocusIndex=t,this.menuFocusIndex>this.menuElement.optionElements.length-1&&(this.menuFocusIndex=this.menuElement.optionElements.length-1),this.menuFocusOptionId=this.menuElement.optionElements[this.menuFocusIndex].id,this.inputElement.setAttribute("aria-owns",this.menuId),this.inputElement.setAttribute("aria-expanded","true"),this.inputElement.setAttribute("aria-activedescendant",this.menuFocusOptionId);const e=this.menuElement.optionElements[this.menuFocusIndex];e.setAttribute("aria-selected","true"),this.menuElement.scrollTo(0,e.offsetTop)}updateListItemTemplate(){var t;this.activeListItemTemplate=null!==(t=this.listItemTemplate)&&void 0!==t?t:this.defaultListItemTemplate}updateOptionTemplate(){var t;this.activeMenuOptionTemplate=null!==(t=this.menuOptionTemplate)&&void 0!==t?t:this.defaultMenuOptionTemplate}updateFilteredOptions(){this.filteredOptionsList=this.optionsList.slice(0),this.filterSelected&&(this.filteredOptionsList=this.filteredOptionsList.filter((t=>-1===this.selectedItems.indexOf(t)))),this.filterQuery&&""!==this.query&&void 0!==this.query&&(this.filteredOptionsList=this.filteredOptionsList.filter((t=>-1!==t.indexOf(this.query))))}updateMenuConfig(){let t=this.configLookup[this.menuPlacement];null===t&&(t=Gt),this.menuConfig=Object.assign(Object.assign({},t),{autoUpdateMode:"auto",fixedPlacement:!0,horizontalViewportLock:!1,verticalViewportLock:!1})}}(0,d.gn)([(0,s.attr)({attribute:"selection"})],gs.prototype,"selection",void 0),(0,d.gn)([(0,s.attr)({attribute:"options"})],gs.prototype,"options",void 0),(0,d.gn)([(0,s.attr)({attribute:"filter-selected",mode:"boolean"})],gs.prototype,"filterSelected",void 0),(0,d.gn)([(0,s.attr)({attribute:"filter-query",mode:"boolean"})],gs.prototype,"filterQuery",void 0),(0,d.gn)([(0,s.attr)({attribute:"max-selected"})],gs.prototype,"maxSelected",void 0),(0,d.gn)([(0,s.attr)({attribute:"no-suggestions-text"})],gs.prototype,"noSuggestionsText",void 0),(0,d.gn)([(0,s.attr)({attribute:"suggestions-available-text"})],gs.prototype,"suggestionsAvailableText",void 0),(0,d.gn)([(0,s.attr)({attribute:"loading-text"})],gs.prototype,"loadingText",void 0),(0,d.gn)([(0,s.attr)({attribute:"label"})],gs.prototype,"label",void 0),(0,d.gn)([(0,s.attr)({attribute:"labelledby"})],gs.prototype,"labelledBy",void 0),(0,d.gn)([(0,s.attr)({attribute:"placeholder"})],gs.prototype,"placeholder",void 0),(0,d.gn)([(0,s.attr)({attribute:"menu-placement"})],gs.prototype,"menuPlacement",void 0),(0,d.gn)([s.observable],gs.prototype,"showLoading",void 0),(0,d.gn)([s.observable],gs.prototype,"listItemTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"defaultListItemTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"activeListItemTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"menuOptionTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"defaultMenuOptionTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"activeMenuOptionTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"listItemContentsTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"menuOptionContentsTemplate",void 0),(0,d.gn)([s.observable],gs.prototype,"optionsList",void 0),(0,d.gn)([s.observable],gs.prototype,"query",void 0),(0,d.gn)([s.observable],gs.prototype,"filteredOptionsList",void 0),(0,d.gn)([s.observable],gs.prototype,"flyoutOpen",void 0),(0,d.gn)([s.observable],gs.prototype,"menuId",void 0),(0,d.gn)([s.observable],gs.prototype,"selectedListTag",void 0),(0,d.gn)([s.observable],gs.prototype,"menuTag",void 0),(0,d.gn)([s.observable],gs.prototype,"menuFocusIndex",void 0),(0,d.gn)([s.observable],gs.prototype,"menuFocusOptionId",void 0),(0,d.gn)([s.observable],gs.prototype,"showNoOptions",void 0),(0,d.gn)([s.observable],gs.prototype,"menuConfig",void 0),(0,d.gn)([s.observable],gs.prototype,"selectedItems",void 0);const bs=(t,e)=>s.html`
        <template role="list" slot="menu-region">
            <div class="options-display" part="options-display">
                <div class="header-region" part="header-region">
                    <slot name="header-region" ${(0,s.slotted)("headerElements")}></slot>
                </div>

                <slot ${(0,s.slotted)("menuElements")}></slot>
                <div class="footer-region" part="footer-region">
                    <slot name="footer-region" ${(0,s.slotted)("footerElements")}></slot>
                </div>
                <div
                    role="alert"
                    aria-live="polite"
                    part="suggestions-available-alert"
                    class="suggestions-available-alert"
                >
                    ${t=>t.suggestionsAvailableText}
                </div>
            </div>
        </template>
    `,fs=(t,e)=>s.html`
        <template
            role="listitem"
            tabindex="-1"
            @click="${(t,e)=>t.handleClick(e.event)}"
        >
            <slot></slot>
        </template>
    `,ys=(t,e)=>s.html`
        <template slot="list-region" role="list" class="picker-list">
            <slot></slot>
            <slot name="input-region"></slot>
        </template>
    `,Cs=(t,e)=>s.html`
        <template
            role="listitem"
            tabindex="0"
            @click="${(t,e)=>t.handleClick(e.event)}"
            @keydown="${(t,e)=>t.handleKeyDown(e.event)}"
        >
            <slot></slot>
        </template>
    `,xs={menuitem:"menuitem",menuitemcheckbox:"menuitemcheckbox",menuitemradio:"menuitemradio"},ws={[xs.menuitem]:"menuitem",[xs.menuitemcheckbox]:"menuitemcheckbox",[xs.menuitemradio]:"menuitemradio"},$s=(t,e)=>s.html`
    <template
        role="${t=>t.role}"
        aria-haspopup="${t=>t.hasSubmenu?"menu":void 0}"
        aria-checked="${t=>t.role!==xs.menuitem?t.checked:void 0}"
        aria-disabled="${t=>t.disabled}"
        aria-expanded="${t=>t.expanded}"
        @keydown="${(t,e)=>t.handleMenuItemKeyDown(e.event)}"
        @click="${(t,e)=>t.handleMenuItemClick(e.event)}"
        @mouseover="${(t,e)=>t.handleMouseOver(e.event)}"
        @mouseout="${(t,e)=>t.handleMouseOut(e.event)}"
        class="${t=>t.disabled?"disabled":""} ${t=>t.expanded?"expanded":""} ${t=>`indent-${t.startColumnCount}`}"
    >
            ${(0,s.when)((t=>t.role===xs.menuitemcheckbox),s.html`
                    <div part="input-container" class="input-container">
                        <span part="checkbox" class="checkbox">
                            <slot name="checkbox-indicator">
                                ${e.checkboxIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
            ${(0,s.when)((t=>t.role===xs.menuitemradio),s.html`
                    <div part="input-container" class="input-container">
                        <span part="radio" class="radio">
                            <slot name="radio-indicator">
                                ${e.radioIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
        </div>
        ${a(t,e)}
        <span class="content" part="content">
            <slot></slot>
        </span>
        ${n(t,e)}
        ${(0,s.when)((t=>t.hasSubmenu),s.html`
                <div
                    part="expand-collapse-glyph-container"
                    class="expand-collapse-glyph-container"
                >
                    <span part="expand-collapse" class="expand-collapse">
                        <slot name="expand-collapse-indicator">
                            ${e.expandCollapseGlyph||""}
                        </slot>
                    </span>
                </div>
            `)}
        ${(0,s.when)((t=>t.expanded),s.html`
                <${t.tagFor(qt)}
                    :anchorElement="${t=>t}"
                    vertical-positioning-mode="dynamic"
                    vertical-default-position="bottom"
                    vertical-inset="true"
                    horizontal-positioning-mode="dynamic"
                    horizontal-default-position="end"
                    class="submenu-region"
                    dir="${t=>t.currentDirection}"
                    @loaded="${t=>t.submenuLoaded()}"
                    ${(0,s.ref)("submenuRegion")}
                    part="submenu-region"
                >
                    <slot name="submenu"></slot>
                </${t.tagFor(qt)}>
            `)}
    </template>
`;class Is extends rt{constructor(){super(...arguments),this.role=xs.menuitem,this.hasSubmenu=!1,this.currentDirection=Lt.N.ltr,this.focusSubmenuOnLoad=!1,this.handleMenuItemKeyDown=t=>{if(t.defaultPrevented)return!1;switch(t.key){case ft:case wt:return this.invoke(),!1;case gt:return this.expandAndFocus(),!1;case vt:if(this.expanded)return this.expanded=!1,this.focus(),!1}return!0},this.handleMenuItemClick=t=>(t.defaultPrevented||this.disabled||this.invoke(),!1),this.submenuLoaded=()=>{this.focusSubmenuOnLoad&&(this.focusSubmenuOnLoad=!1,this.hasSubmenu&&(this.submenu.focus(),this.setAttribute("tabindex","-1")))},this.handleMouseOver=t=>(this.disabled||!this.hasSubmenu||this.expanded||(this.expanded=!0),!1),this.handleMouseOut=t=>(!this.expanded||this.contains(document.activeElement)||(this.expanded=!1),!1),this.expandAndFocus=()=>{this.hasSubmenu&&(this.focusSubmenuOnLoad=!0,this.expanded=!0)},this.invoke=()=>{if(!this.disabled)switch(this.role){case xs.menuitemcheckbox:this.checked=!this.checked;break;case xs.menuitem:this.updateSubmenu(),this.hasSubmenu?this.expandAndFocus():this.$emit("change");break;case xs.menuitemradio:this.checked||(this.checked=!0)}},this.updateSubmenu=()=>{this.submenu=this.domChildren().find((t=>"menu"===t.getAttribute("role"))),this.hasSubmenu=void 0!==this.submenu}}expandedChanged(t){if(this.$fastController.isConnected){if(void 0===this.submenu)return;!1===this.expanded?this.submenu.collapseExpandedItem():this.currentDirection=Bt(this),this.$emit("expanded-change",this,{bubbles:!1})}}checkedChanged(t,e){this.$fastController.isConnected&&this.$emit("change")}connectedCallback(){super.connectedCallback(),s.DOM.queueUpdate((()=>{this.updateSubmenu()})),this.startColumnCount||(this.startColumnCount=1),this.observer=new MutationObserver(this.updateSubmenu)}disconnectedCallback(){super.disconnectedCallback(),this.submenu=void 0,void 0!==this.observer&&(this.observer.disconnect(),this.observer=void 0)}domChildren(){return Array.from(this.children).filter((t=>!t.hasAttribute("hidden")))}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],Is.prototype,"disabled",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Is.prototype,"expanded",void 0),(0,d.gn)([s.observable],Is.prototype,"startColumnCount",void 0),(0,d.gn)([s.attr],Is.prototype,"role",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Is.prototype,"checked",void 0),(0,d.gn)([s.observable],Is.prototype,"submenuRegion",void 0),(0,d.gn)([s.observable],Is.prototype,"hasSubmenu",void 0),(0,d.gn)([s.observable],Is.prototype,"currentDirection",void 0),(0,d.gn)([s.observable],Is.prototype,"submenu",void 0),dt(Is,o);const ks=(t,e)=>s.html`
    <template
        slot="${t=>t.slot?t.slot:t.isNestedMenu()?"submenu":void 0}"
        role="menu"
        @keydown="${(t,e)=>t.handleMenuKeyDown(e.event)}"
        @focusout="${(t,e)=>t.handleFocusOut(e.event)}"
    >
        <slot ${(0,s.slotted)("items")}></slot>
    </template>
`;class Es extends rt{constructor(){super(...arguments),this.expandedItem=null,this.focusIndex=-1,this.isNestedMenu=()=>null!==this.parentElement&&Be(this.parentElement)&&"menuitem"===this.parentElement.getAttribute("role"),this.handleFocusOut=t=>{if(!this.contains(t.relatedTarget)&&void 0!==this.menuItems){this.collapseExpandedItem();const t=this.menuItems.findIndex(this.isFocusableElement);this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.menuItems[t].setAttribute("tabindex","0"),this.focusIndex=t}},this.handleItemFocus=t=>{const e=t.target;void 0!==this.menuItems&&e!==this.menuItems[this.focusIndex]&&(this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.focusIndex=this.menuItems.indexOf(e),e.setAttribute("tabindex","0"))},this.handleExpandedChanged=t=>{if(t.defaultPrevented||null===t.target||void 0===this.menuItems||this.menuItems.indexOf(t.target)<0)return;t.preventDefault();const e=t.target;null===this.expandedItem||e!==this.expandedItem||!1!==e.expanded?e.expanded&&(null!==this.expandedItem&&this.expandedItem!==e&&(this.expandedItem.expanded=!1),this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.expandedItem=e,this.focusIndex=this.menuItems.indexOf(e),e.setAttribute("tabindex","0")):this.expandedItem=null},this.removeItemListeners=()=>{void 0!==this.menuItems&&this.menuItems.forEach((t=>{t.removeEventListener("expanded-change",this.handleExpandedChanged),t.removeEventListener("focus",this.handleItemFocus)}))},this.setItems=()=>{const t=this.domChildren();this.removeItemListeners(),this.menuItems=t;const e=this.menuItems.filter(this.isMenuItemElement);e.length&&(this.focusIndex=0);const i=e.reduce(((t,e)=>{const i=function(t){const e=t.getAttribute("role"),i=t.querySelector("[slot=start]");return e!==xs.menuitem&&null===i||e===xs.menuitem&&null!==i?1:e!==xs.menuitem&&null!==i?2:0}(e);return t>i?t:i}),0);e.forEach(((t,e)=>{t.setAttribute("tabindex",0===e?"0":"-1"),t.addEventListener("expanded-change",this.handleExpandedChanged),t.addEventListener("focus",this.handleItemFocus),t instanceof Is&&(t.startColumnCount=i)}))},this.changeHandler=t=>{if(void 0===this.menuItems)return;const e=t.target,i=this.menuItems.indexOf(e);if(-1!==i&&"menuitemradio"===e.role&&!0===e.checked){for(let t=i-1;t>=0;--t){const e=this.menuItems[t],i=e.getAttribute("role");if(i===xs.menuitemradio&&(e.checked=!1),"separator"===i)break}const t=this.menuItems.length-1;for(let e=i+1;e<=t;++e){const t=this.menuItems[e],i=t.getAttribute("role");if(i===xs.menuitemradio&&(t.checked=!1),"separator"===i)break}}},this.isMenuItemElement=t=>Be(t)&&Es.focusableElementRoles.hasOwnProperty(t.getAttribute("role")),this.isFocusableElement=t=>this.isMenuItemElement(t)}itemsChanged(t,e){this.$fastController.isConnected&&void 0!==this.menuItems&&this.setItems()}connectedCallback(){super.connectedCallback(),s.DOM.queueUpdate((()=>{this.setItems()})),this.addEventListener("change",this.changeHandler)}disconnectedCallback(){super.disconnectedCallback(),this.removeItemListeners(),this.menuItems=void 0,this.removeEventListener("change",this.changeHandler)}focus(){this.setFocus(0,1)}collapseExpandedItem(){null!==this.expandedItem&&(this.expandedItem.expanded=!1,this.expandedItem=null)}handleMenuKeyDown(t){if(!t.defaultPrevented&&void 0!==this.menuItems)switch(t.key){case mt:return void this.setFocus(this.focusIndex+1,1);case bt:return void this.setFocus(this.focusIndex-1,-1);case xt:return void this.setFocus(this.menuItems.length-1,-1);case Ct:return void this.setFocus(0,1);default:return!0}}domChildren(){return Array.from(this.children).filter((t=>!t.hasAttribute("hidden")))}setFocus(t,e){if(void 0!==this.menuItems)for(;t>=0&&t<this.menuItems.length;){const i=this.menuItems[t];if(this.isFocusableElement(i)){this.focusIndex>-1&&this.menuItems.length>=this.focusIndex-1&&this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.focusIndex=t,i.setAttribute("tabindex","0"),i.focus();break}t+=e}}}Es.focusableElementRoles=ws,(0,d.gn)([s.observable],Es.prototype,"items",void 0);const Ts=(t,e)=>s.html`
    <template class="${t=>t.readOnly?"readonly":""}">
        <label
            part="label"
            for="control"
            class="${t=>t.defaultSlottedNodes&&t.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
        <div class="root" part="root">
            ${a(t,e)}
            <input
                class="control"
                part="control"
                id="control"
                @input="${t=>t.handleTextInput()}"
                @change="${t=>t.handleChange()}"
                @keydown="${(t,e)=>t.handleKeyDown(e.event)}"
                @blur="${(t,e)=>t.handleBlur()}"
                ?autofocus="${t=>t.autofocus}"
                ?disabled="${t=>t.disabled}"
                list="${t=>t.list}"
                maxlength="${t=>t.maxlength}"
                minlength="${t=>t.minlength}"
                placeholder="${t=>t.placeholder}"
                ?readonly="${t=>t.readOnly}"
                ?required="${t=>t.required}"
                size="${t=>t.size}"
                type="text"
                inputmode="numeric"
                min="${t=>t.min}"
                max="${t=>t.max}"
                step="${t=>t.step}"
                aria-atomic="${t=>t.ariaAtomic}"
                aria-busy="${t=>t.ariaBusy}"
                aria-controls="${t=>t.ariaControls}"
                aria-current="${t=>t.ariaCurrent}"
                aria-describedby="${t=>t.ariaDescribedby}"
                aria-details="${t=>t.ariaDetails}"
                aria-disabled="${t=>t.ariaDisabled}"
                aria-errormessage="${t=>t.ariaErrormessage}"
                aria-flowto="${t=>t.ariaFlowto}"
                aria-haspopup="${t=>t.ariaHaspopup}"
                aria-hidden="${t=>t.ariaHidden}"
                aria-invalid="${t=>t.ariaInvalid}"
                aria-keyshortcuts="${t=>t.ariaKeyshortcuts}"
                aria-label="${t=>t.ariaLabel}"
                aria-labelledby="${t=>t.ariaLabelledby}"
                aria-live="${t=>t.ariaLive}"
                aria-owns="${t=>t.ariaOwns}"
                aria-relevant="${t=>t.ariaRelevant}"
                aria-roledescription="${t=>t.ariaRoledescription}"
                ${(0,s.ref)("control")}
            />
            ${(0,s.when)((t=>!t.hideStep&&!t.readOnly&&!t.disabled),s.html`
                    <div class="controls" part="controls">
                        <div class="step-up" part="step-up" @click="${t=>t.stepUp()}">
                            <slot name="step-up-glyph">
                                ${e.stepUpGlyph||""}
                            </slot>
                        </div>
                        <div
                            class="step-down"
                            part="step-down"
                            @click="${t=>t.stepDown()}"
                        >
                            <slot name="step-down-glyph">
                                ${e.stepDownGlyph||""}
                            </slot>
                        </div>
                    </div>
                `)}
            ${n(t,e)}
        </div>
    </template>
`;class Os extends rt{}class Ds extends(he(Os)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}const Rs={email:"email",password:"password",tel:"tel",text:"text",url:"url"};class Ss extends Ds{constructor(){super(...arguments),this.type=Rs.text}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}typeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.type=this.type,this.validate())}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.maxLength=this.maxlength,this.validate())}minlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.minLength=this.minlength,this.validate())}patternChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.pattern=this.pattern,this.validate())}sizeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.size=this.size)}spellcheckChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.spellcheck=this.spellcheck)}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type",this.type),this.validate(),this.autofocus&&s.DOM.queueUpdate((()=>{this.focus()}))}select(){this.control.select(),this.$emit("select")}handleTextInput(){this.value=this.control.value}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Ss.prototype,"readOnly",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Ss.prototype,"autofocus",void 0),(0,d.gn)([s.attr],Ss.prototype,"placeholder",void 0),(0,d.gn)([s.attr],Ss.prototype,"type",void 0),(0,d.gn)([s.attr],Ss.prototype,"list",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ss.prototype,"maxlength",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ss.prototype,"minlength",void 0),(0,d.gn)([s.attr],Ss.prototype,"pattern",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ss.prototype,"size",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Ss.prototype,"spellcheck",void 0),(0,d.gn)([s.observable],Ss.prototype,"defaultSlottedNodes",void 0);class As{}dt(As,Rt),dt(Ss,o,As);class Fs extends rt{}class Ls extends(he(Fs)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class Ps extends Ls{constructor(){super(...arguments),this.hideStep=!1,this.step=1,this.isUserInput=!1}maxChanged(t,e){var i;this.max=Math.max(e,null!==(i=this.min)&&void 0!==i?i:e);const s=Math.min(this.min,this.max);void 0!==this.min&&this.min!==s&&(this.min=s),this.value=this.getValidValue(this.value)}minChanged(t,e){var i;this.min=Math.min(e,null!==(i=this.max)&&void 0!==i?i:e);const s=Math.max(this.min,this.max);void 0!==this.max&&this.max!==s&&(this.max=s),this.value=this.getValidValue(this.value)}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(t){this.value=t.toString()}valueChanged(t,e){this.value=this.getValidValue(e),e===this.value&&(this.control&&!this.isUserInput&&(this.control.value=this.value),super.valueChanged(t,this.value),void 0===t||this.isUserInput||(this.$emit("input"),this.$emit("change")),this.isUserInput=!1)}validate(){super.validate(this.control)}getValidValue(t){var e,i;let s=parseFloat(parseFloat(t).toPrecision(12));return isNaN(s)?s="":(s=Math.min(s,null!==(e=this.max)&&void 0!==e?e:s),s=Math.max(s,null!==(i=this.min)&&void 0!==i?i:s).toString()),s}stepUp(){const t=parseFloat(this.value),e=isNaN(t)?this.min>0?this.min:this.max<0?this.max:this.min?0:this.step:t+this.step;this.value=e.toString()}stepDown(){const t=parseFloat(this.value),e=isNaN(t)?this.min>0?this.min:this.max<0?this.max:this.min?0:0-this.step:t-this.step;this.value=e.toString()}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type","number"),this.validate(),this.control.value=this.value,this.autofocus&&s.DOM.queueUpdate((()=>{this.focus()}))}select(){this.control.select(),this.$emit("select")}handleTextInput(){this.control.value=this.control.value.replace(/[^0-9\-+e.]/g,""),this.isUserInput=!0,this.value=this.control.value}handleChange(){this.$emit("change")}handleKeyDown(t){switch(t.key){case bt:return this.stepUp(),!1;case mt:return this.stepDown(),!1}return!0}handleBlur(){this.control.value=this.value}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Ps.prototype,"readOnly",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Ps.prototype,"autofocus",void 0),(0,d.gn)([(0,s.attr)({attribute:"hide-step",mode:"boolean"})],Ps.prototype,"hideStep",void 0),(0,d.gn)([s.attr],Ps.prototype,"placeholder",void 0),(0,d.gn)([s.attr],Ps.prototype,"list",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ps.prototype,"maxlength",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ps.prototype,"minlength",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ps.prototype,"size",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ps.prototype,"step",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ps.prototype,"max",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Ps.prototype,"min",void 0),(0,d.gn)([s.observable],Ps.prototype,"defaultSlottedNodes",void 0),dt(Ps,o,As);const Ms=(t,e)=>s.html`
    <template
        role="progressbar"
        aria-valuenow="${t=>t.value}"
        aria-valuemin="${t=>t.min}"
        aria-valuemax="${t=>t.max}"
        class="${t=>t.paused?"paused":""}"
    >
        ${(0,s.when)((t=>"number"==typeof t.value),s.html`
                <svg
                    class="progress"
                    part="progress"
                    viewBox="0 0 16 16"
                    slot="determinate"
                >
                    <circle
                        class="background"
                        part="background"
                        cx="8px"
                        cy="8px"
                        r="7px"
                    ></circle>
                    <circle
                        class="determinate"
                        part="determinate"
                        style="stroke-dasharray: ${t=>44*t.percentComplete/100}px ${44}px"
                        cx="8px"
                        cy="8px"
                        r="7px"
                    ></circle>
                </svg>
            `,s.html`
                <slot name="indeterminate" slot="indeterminate">
                    ${e.indeterminateIndicator||""}
                </slot>
            `)}
    </template>
`;class Hs extends rt{constructor(){super(...arguments),this.percentComplete=0}valueChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}minChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}maxChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}connectedCallback(){super.connectedCallback(),this.updatePercentComplete()}updatePercentComplete(){const t="number"==typeof this.min?this.min:0,e="number"==typeof this.max?this.max:100,i="number"==typeof this.value?this.value:0,s=e-t;this.percentComplete=0===s?0:Math.fround((i-t)/s*100)}}(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Hs.prototype,"value",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Hs.prototype,"min",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Hs.prototype,"max",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Hs.prototype,"paused",void 0),(0,d.gn)([s.observable],Hs.prototype,"percentComplete",void 0);const zs=(t,e)=>s.html`
    <template
        role="progressbar"
        aria-valuenow="${t=>t.value}"
        aria-valuemin="${t=>t.min}"
        aria-valuemax="${t=>t.max}"
        class="${t=>t.paused?"paused":""}"
    >
        ${(0,s.when)((t=>"number"==typeof t.value),s.html`
                <div class="progress" part="progress" slot="determinate">
                    <div
                        class="determinate"
                        part="determinate"
                        style="width: ${t=>t.percentComplete}%"
                    ></div>
                </div>
            `,s.html`
                <div class="progress" part="progress" slot="indeterminate">
                    <slot class="indeterminate" name="indeterminate">
                        ${e.indeterminateIndicator1||""}
                        ${e.indeterminateIndicator2||""}
                    </slot>
                </div>
            `)}
    </template>
`,Vs=(t,e)=>s.html`
    <template
        role="radiogroup"
        aria-disabled="${t=>t.disabled}"
        aria-readonly="${t=>t.readOnly}"
        @click="${(t,e)=>t.clickHandler(e.event)}"
        @keydown="${(t,e)=>t.keydownHandler(e.event)}"
        @focusout="${(t,e)=>t.focusOutHandler(e.event)}"
    >
        <slot name="label"></slot>
        <div
            class="positioning-region ${t=>t.orientation===Qi.i.horizontal?"horizontal":"vertical"}"
            part="positioning-region"
        >
            <slot
                ${(0,s.slotted)({property:"slottedRadioButtons",filter:(0,s.elements)("[role=radio]")})}
            ></slot>
        </div>
    </template>
`;class Ns extends rt{constructor(){super(...arguments),this.orientation=Qi.i.horizontal,this.radioChangeHandler=t=>{const e=t.target;e.checked&&(this.slottedRadioButtons.forEach((t=>{t!==e&&(t.checked=!1,this.isInsideFoundationToolbar||t.setAttribute("tabindex","-1"))})),this.selectedRadio=e,this.value=e.value,e.setAttribute("tabindex","0"),this.focusedRadio=e),t.stopPropagation()},this.moveToRadioByIndex=(t,e)=>{const i=t[e];this.isInsideToolbar||(i.setAttribute("tabindex","0"),i.readOnly?this.slottedRadioButtons.forEach((t=>{t!==i&&t.setAttribute("tabindex","-1")})):(i.checked=!0,this.selectedRadio=i)),this.focusedRadio=i,i.focus()},this.moveRightOffGroup=()=>{var t;null===(t=this.nextElementSibling)||void 0===t||t.focus()},this.moveLeftOffGroup=()=>{var t;null===(t=this.previousElementSibling)||void 0===t||t.focus()},this.focusOutHandler=t=>{const e=this.slottedRadioButtons,i=t.target,s=null!==i?e.indexOf(i):0,o=this.focusedRadio?e.indexOf(this.focusedRadio):-1;return(0===o&&s===o||o===e.length-1&&o===s)&&(this.selectedRadio?(this.focusedRadio=this.selectedRadio,this.isInsideFoundationToolbar||(this.selectedRadio.setAttribute("tabindex","0"),e.forEach((t=>{t!==this.selectedRadio&&t.setAttribute("tabindex","-1")})))):(this.focusedRadio=e[0],this.focusedRadio.setAttribute("tabindex","0"),e.forEach((t=>{t!==this.focusedRadio&&t.setAttribute("tabindex","-1")})))),!0},this.clickHandler=t=>{const e=t.target;if(e){const t=this.slottedRadioButtons;e.checked||0===t.indexOf(e)?(e.setAttribute("tabindex","0"),this.selectedRadio=e):(e.setAttribute("tabindex","-1"),this.selectedRadio=null),this.focusedRadio=e}t.preventDefault()},this.shouldMoveOffGroupToTheRight=(t,e,i)=>t===e.length&&this.isInsideToolbar&&i===gt,this.shouldMoveOffGroupToTheLeft=(t,e)=>(this.focusedRadio?t.indexOf(this.focusedRadio)-1:0)<0&&this.isInsideToolbar&&e===vt,this.checkFocusedRadio=()=>{null===this.focusedRadio||this.focusedRadio.readOnly||this.focusedRadio.checked||(this.focusedRadio.checked=!0,this.focusedRadio.setAttribute("tabindex","0"),this.focusedRadio.focus(),this.selectedRadio=this.focusedRadio)},this.moveRight=t=>{const e=this.slottedRadioButtons;let i=0;if(i=this.focusedRadio?e.indexOf(this.focusedRadio)+1:1,this.shouldMoveOffGroupToTheRight(i,e,t.key))this.moveRightOffGroup();else for(i===e.length&&(i=0);i<e.length&&e.length>1;){if(!e[i].disabled){this.moveToRadioByIndex(e,i);break}if(this.focusedRadio&&i===e.indexOf(this.focusedRadio))break;if(i+1>=e.length){if(this.isInsideToolbar)break;i=0}else i+=1}},this.moveLeft=t=>{const e=this.slottedRadioButtons;let i=0;if(i=this.focusedRadio?e.indexOf(this.focusedRadio)-1:0,i=i<0?e.length-1:i,this.shouldMoveOffGroupToTheLeft(e,t.key))this.moveLeftOffGroup();else for(;i>=0&&e.length>1;){if(!e[i].disabled){this.moveToRadioByIndex(e,i);break}if(this.focusedRadio&&i===e.indexOf(this.focusedRadio))break;i-1<0?i=e.length-1:i-=1}},this.keydownHandler=t=>{const e=t.key;if(e in It&&this.isInsideFoundationToolbar)return!0;switch(e){case ft:this.checkFocusedRadio();break;case gt:case mt:this.direction===Lt.N.ltr?this.moveRight(t):this.moveLeft(t);break;case vt:case bt:this.direction===Lt.N.ltr?this.moveLeft(t):this.moveRight(t);break;default:return!0}}}readOnlyChanged(){void 0!==this.slottedRadioButtons&&this.slottedRadioButtons.forEach((t=>{this.readOnly?t.readOnly=!0:t.readOnly=!1}))}disabledChanged(){void 0!==this.slottedRadioButtons&&this.slottedRadioButtons.forEach((t=>{this.disabled?t.disabled=!0:t.disabled=!1}))}nameChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach((t=>{t.setAttribute("name",this.name)}))}valueChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach((t=>{t.value===this.value&&(t.checked=!0,this.selectedRadio=t)})),this.$emit("change")}slottedRadioButtonsChanged(t,e){this.slottedRadioButtons&&this.slottedRadioButtons.length>0&&this.setupRadioButtons()}get parentToolbar(){return this.closest('[role="toolbar"]')}get isInsideToolbar(){var t;return null!==(t=this.parentToolbar)&&void 0!==t&&t}get isInsideFoundationToolbar(){var t;return!!(null===(t=this.parentToolbar)||void 0===t?void 0:t.$fastController)}connectedCallback(){super.connectedCallback(),this.direction=Bt(this),this.setupRadioButtons()}disconnectedCallback(){this.slottedRadioButtons.forEach((t=>{t.removeEventListener("change",this.radioChangeHandler)}))}setupRadioButtons(){const t=this.slottedRadioButtons.filter((t=>t.hasAttribute("checked"))),e=t?t.length:0;e>1&&(t[e-1].checked=!0);let i=!1;if(this.slottedRadioButtons.forEach((t=>{void 0!==this.name&&t.setAttribute("name",this.name),this.disabled&&(t.disabled=!0),this.readOnly&&(t.readOnly=!0),this.value&&this.value===t.value?(this.selectedRadio=t,this.focusedRadio=t,t.checked=!0,t.setAttribute("tabindex","0"),i=!0):(this.isInsideFoundationToolbar||t.setAttribute("tabindex","-1"),t.checked=!1),t.addEventListener("change",this.radioChangeHandler)})),void 0===this.value&&this.slottedRadioButtons.length>0){const t=this.slottedRadioButtons.filter((t=>t.hasAttribute("checked"))),e=null!==t?t.length:0;if(e>0&&!i){const i=t[e-1];i.checked=!0,this.focusedRadio=i,i.setAttribute("tabindex","0")}else this.slottedRadioButtons[0].setAttribute("tabindex","0"),this.focusedRadio=this.slottedRadioButtons[0]}}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Ns.prototype,"readOnly",void 0),(0,d.gn)([(0,s.attr)({attribute:"disabled",mode:"boolean"})],Ns.prototype,"disabled",void 0),(0,d.gn)([s.attr],Ns.prototype,"name",void 0),(0,d.gn)([s.attr],Ns.prototype,"value",void 0),(0,d.gn)([s.attr],Ns.prototype,"orientation",void 0),(0,d.gn)([s.observable],Ns.prototype,"childItems",void 0),(0,d.gn)([s.observable],Ns.prototype,"slottedRadioButtons",void 0);const Bs=(t,e)=>s.html`
    <template
        role="radio"
        class="${t=>t.checked?"checked":""} ${t=>t.readOnly?"readonly":""}"
        aria-checked="${t=>t.checked}"
        aria-required="${t=>t.required}"
        aria-disabled="${t=>t.disabled}"
        aria-readonly="${t=>t.readOnly}"
        @keypress="${(t,e)=>t.keypressHandler(e.event)}"
        @click="${(t,e)=>t.clickHandler(e.event)}"
    >
        <div part="control" class="control">
            <slot name="checked-indicator">
                ${e.checkedIndicator||""}
            </slot>
        </div>
        <label
            part="label"
            class="${t=>t.defaultSlottedNodes&&t.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
    </template>
`;class qs extends rt{}class Us extends(de(qs)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class js extends Us{constructor(){super(),this.initialValue="on",this.keypressHandler=t=>{if(t.key!==wt)return!0;this.checked||this.readOnly||(this.checked=!0)},this.proxy.setAttribute("type","radio")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}defaultCheckedChanged(){var t;this.$fastController.isConnected&&!this.dirtyChecked&&(this.isInsideRadioGroup()||(this.checked=null!==(t=this.defaultChecked)&&void 0!==t&&t,this.dirtyChecked=!1))}connectedCallback(){var t,e;super.connectedCallback(),this.validate(),"radiogroup"!==(null===(t=this.parentElement)||void 0===t?void 0:t.getAttribute("role"))&&null===this.getAttribute("tabindex")&&(this.disabled||this.setAttribute("tabindex","0")),this.checkedAttribute&&(this.dirtyChecked||this.isInsideRadioGroup()||(this.checked=null!==(e=this.defaultChecked)&&void 0!==e&&e,this.dirtyChecked=!1))}isInsideRadioGroup(){return null!==this.closest("[role=radiogroup]")}clickHandler(t){this.disabled||this.readOnly||this.checked||(this.checked=!0)}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],js.prototype,"readOnly",void 0),(0,d.gn)([s.observable],js.prototype,"name",void 0),(0,d.gn)([s.observable],js.prototype,"defaultSlottedNodes",void 0);class _s extends rt{constructor(){super(...arguments),this.framesPerSecond=60,this.updatingItems=!1,this.speed=600,this.easing="ease-in-out",this.flippersHiddenFromAT=!1,this.scrolling=!1,this.resizeDetector=null}get frameTime(){return 1e3/this.framesPerSecond}scrollingChanged(t,e){if(this.scrollContainer){const t=1==this.scrolling?"scrollstart":"scrollend";this.$emit(t,this.scrollContainer.scrollLeft)}}get isRtl(){return this.scrollItems.length>1&&this.scrollItems[0].offsetLeft>this.scrollItems[1].offsetLeft}connectedCallback(){super.connectedCallback(),this.initializeResizeDetector()}disconnectedCallback(){this.disconnectResizeDetector(),super.disconnectedCallback()}scrollItemsChanged(t,e){e&&!this.updatingItems&&s.DOM.queueUpdate((()=>this.setStops()))}disconnectResizeDetector(){this.resizeDetector&&(this.resizeDetector.disconnect(),this.resizeDetector=null)}initializeResizeDetector(){this.disconnectResizeDetector(),this.resizeDetector=new window.ResizeObserver(this.resized.bind(this)),this.resizeDetector.observe(this)}updateScrollStops(){this.updatingItems=!0;const t=this.scrollItems.reduce(((t,e)=>e instanceof HTMLSlotElement?t.concat(e.assignedElements()):(t.push(e),t)),[]);this.scrollItems=t,this.updatingItems=!1}setStops(){this.updateScrollStops();const{scrollContainer:t}=this,{scrollLeft:e}=t,{width:i,left:s}=t.getBoundingClientRect();this.width=i;let o=0,n=this.scrollItems.map(((t,i)=>{const{left:n,width:a}=t.getBoundingClientRect(),r=Math.round(n+e-s),l=Math.round(r+a);return this.isRtl?-l:(o=l,0===i?0:r)})).concat(o);n=this.fixScrollMisalign(n),n.sort(((t,e)=>Math.abs(t)-Math.abs(e))),this.scrollStops=n,this.setFlippers()}validateStops(t=!0){const e=()=>!!this.scrollStops.find((t=>t>0));return!e()&&t&&this.setStops(),e()}fixScrollMisalign(t){if(this.isRtl&&t.some((t=>t>0))){t.sort(((t,e)=>e-t));const e=t[0];t=t.map((t=>t-e))}return t}setFlippers(){var t,e;const i=this.scrollContainer.scrollLeft;if(null===(t=this.previousFlipperContainer)||void 0===t||t.classList.toggle("disabled",0===i),this.scrollStops){const t=Math.abs(this.scrollStops[this.scrollStops.length-1]);null===(e=this.nextFlipperContainer)||void 0===e||e.classList.toggle("disabled",this.validateStops(!1)&&Math.abs(i)+this.width>=t)}}scrollInView(t,e=0,i){var s;if("number"!=typeof t&&t&&(t=this.scrollItems.findIndex((e=>e===t||e.contains(t)))),void 0!==t){i=null!=i?i:e;const{scrollContainer:o,scrollStops:n,scrollItems:a}=this,{scrollLeft:r}=this.scrollContainer,{width:l}=o.getBoundingClientRect(),h=n[t],{width:d}=a[t].getBoundingClientRect(),c=h+d,u=r+e>h;if(u||r+l-i<c){const t=null!==(s=[...n].sort(((t,e)=>u?e-t:t-e)).find((t=>u?t+e<h:t+l-(null!=i?i:0)>c)))&&void 0!==s?s:0;this.scrollToPosition(t)}}}keyupHandler(t){switch(t.key){case"ArrowLeft":this.scrollToPrevious();break;case"ArrowRight":this.scrollToNext()}}scrollToPrevious(){this.validateStops();const t=this.scrollContainer.scrollLeft,e=this.scrollStops.findIndex(((e,i)=>e>=t&&(this.isRtl||i===this.scrollStops.length-1||this.scrollStops[i+1]>t))),i=Math.abs(this.scrollStops[e+1]);let s=this.scrollStops.findIndex((t=>Math.abs(t)+this.width>i));(s>=e||-1===s)&&(s=e>0?e-1:0),this.scrollToPosition(this.scrollStops[s],t)}scrollToNext(){this.validateStops();const t=this.scrollContainer.scrollLeft,e=this.scrollStops.findIndex((e=>Math.abs(e)>=Math.abs(t))),i=this.scrollStops.findIndex((e=>Math.abs(t)+this.width<=Math.abs(e)));let s=e;i>e+2?s=i-2:e<this.scrollStops.length-2&&(s=e+1),this.scrollToPosition(this.scrollStops[s],t)}scrollToPosition(t,e=this.scrollContainer.scrollLeft){var i;if(this.scrolling)return;this.scrolling=!0;const s=null!==(i=this.duration)&&void 0!==i?i:Math.abs(t-e)/this.speed+"s";this.content.style.setProperty("transition-duration",s);const o=parseFloat(getComputedStyle(this.content).getPropertyValue("transition-duration")),n=e=>{e&&e.target!==e.currentTarget||(this.content.style.setProperty("transition-duration","0s"),this.content.style.removeProperty("transform"),this.scrollContainer.style.setProperty("scroll-behavior","auto"),this.scrollContainer.scrollLeft=t,this.setFlippers(),this.content.removeEventListener("transitionend",n),this.scrolling=!1)};if(0===o)return void n();this.content.addEventListener("transitionend",n);const a=this.scrollContainer.scrollWidth-this.scrollContainer.clientWidth;let r=this.scrollContainer.scrollLeft-Math.min(t,a);this.isRtl&&(r=this.scrollContainer.scrollLeft+Math.min(Math.abs(t),a)),this.content.style.setProperty("transition-property","transform"),this.content.style.setProperty("transition-timing-function",this.easing),this.content.style.setProperty("transform",`translateX(${r}px)`)}resized(){this.resizeTimeout&&(this.resizeTimeout=clearTimeout(this.resizeTimeout)),this.resizeTimeout=setTimeout((()=>{this.width=this.scrollContainer.offsetWidth,this.setFlippers()}),this.frameTime)}scrolled(){this.scrollTimeout&&(this.scrollTimeout=clearTimeout(this.scrollTimeout)),this.scrollTimeout=setTimeout((()=>{this.setFlippers()}),this.frameTime)}}(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],_s.prototype,"speed",void 0),(0,d.gn)([s.attr],_s.prototype,"duration",void 0),(0,d.gn)([s.attr],_s.prototype,"easing",void 0),(0,d.gn)([(0,s.attr)({attribute:"flippers-hidden-from-at",converter:s.booleanConverter})],_s.prototype,"flippersHiddenFromAT",void 0),(0,d.gn)([s.observable],_s.prototype,"scrolling",void 0),(0,d.gn)([s.observable],_s.prototype,"scrollItems",void 0),(0,d.gn)([(0,s.attr)({attribute:"view"})],_s.prototype,"view",void 0);const Ks=(t,e)=>{var i,o;return s.html`
    <template
        class="horizontal-scroll"
        @keyup="${(t,e)=>t.keyupHandler(e.event)}"
    >
        ${a(t,e)}
        <div class="scroll-area" part="scroll-area">
            <div
                class="scroll-view"
                part="scroll-view"
                @scroll="${t=>t.scrolled()}"
                ${(0,s.ref)("scrollContainer")}
            >
                <div class="content-container" part="content-container" ${(0,s.ref)("content")}>
                    <slot
                        ${(0,s.slotted)({property:"scrollItems",filter:(0,s.elements)()})}
                    ></slot>
                </div>
            </div>
            ${(0,s.when)((t=>"mobile"!==t.view),s.html`
                    <div
                        class="scroll scroll-prev"
                        part="scroll-prev"
                        ${(0,s.ref)("previousFlipperContainer")}
                    >
                        <div class="scroll-action" part="scroll-action-previous">
                            <slot name="previous-flipper">
                                ${e.previousFlipper instanceof Function?e.previousFlipper(t,e):null!==(i=e.previousFlipper)&&void 0!==i?i:""}
                            </slot>
                        </div>
                    </div>
                    <div
                        class="scroll scroll-next"
                        part="scroll-next"
                        ${(0,s.ref)("nextFlipperContainer")}
                    >
                        <div class="scroll-action" part="scroll-action-next">
                            <slot name="next-flipper">
                                ${e.nextFlipper instanceof Function?e.nextFlipper(t,e):null!==(o=e.nextFlipper)&&void 0!==o?o:""}
                            </slot>
                        </div>
                    </div>
                `)}
        </div>
        ${n(t,e)}
    </template>
`};function Ws(t,e,i){return t.nodeType!==Node.TEXT_NODE||"string"==typeof t.nodeValue&&!!t.nodeValue.trim().length}const Gs=(t,e)=>s.html`
    <template
        class="
            ${t=>t.readOnly?"readonly":""}
        "
    >
        <label
            part="label"
            for="control"
            class="${t=>t.defaultSlottedNodes&&t.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot
                ${(0,s.slotted)({property:"defaultSlottedNodes",filter:Ws})}
            ></slot>
        </label>
        <div class="root" part="root" ${(0,s.ref)("root")}>
            ${a(t,e)}
            <div class="input-wrapper" part="input-wrapper">
                <input
                    class="control"
                    part="control"
                    id="control"
                    @input="${t=>t.handleTextInput()}"
                    @change="${t=>t.handleChange()}"
                    ?autofocus="${t=>t.autofocus}"
                    ?disabled="${t=>t.disabled}"
                    list="${t=>t.list}"
                    maxlength="${t=>t.maxlength}"
                    minlength="${t=>t.minlength}"
                    pattern="${t=>t.pattern}"
                    placeholder="${t=>t.placeholder}"
                    ?readonly="${t=>t.readOnly}"
                    ?required="${t=>t.required}"
                    size="${t=>t.size}"
                    ?spellcheck="${t=>t.spellcheck}"
                    :value="${t=>t.value}"
                    type="search"
                    aria-atomic="${t=>t.ariaAtomic}"
                    aria-busy="${t=>t.ariaBusy}"
                    aria-controls="${t=>t.ariaControls}"
                    aria-current="${t=>t.ariaCurrent}"
                    aria-describedby="${t=>t.ariaDescribedby}"
                    aria-details="${t=>t.ariaDetails}"
                    aria-disabled="${t=>t.ariaDisabled}"
                    aria-errormessage="${t=>t.ariaErrormessage}"
                    aria-flowto="${t=>t.ariaFlowto}"
                    aria-haspopup="${t=>t.ariaHaspopup}"
                    aria-hidden="${t=>t.ariaHidden}"
                    aria-invalid="${t=>t.ariaInvalid}"
                    aria-keyshortcuts="${t=>t.ariaKeyshortcuts}"
                    aria-label="${t=>t.ariaLabel}"
                    aria-labelledby="${t=>t.ariaLabelledby}"
                    aria-live="${t=>t.ariaLive}"
                    aria-owns="${t=>t.ariaOwns}"
                    aria-relevant="${t=>t.ariaRelevant}"
                    aria-roledescription="${t=>t.ariaRoledescription}"
                    ${(0,s.ref)("control")}
                />
                <slot name="close-button">
                    <button
                        class="clear-button ${t=>t.value?"":"clear-button__hidden"}"
                        part="clear-button"
                        tabindex="-1"
                        @click=${t=>t.handleClearInput()}
                    >
                        <slot name="close-glyph">
                            <svg
                                width="9"
                                height="9"
                                viewBox="0 0 9 9"
                                xmlns="http://www.w3.org/2000/svg"
                            >
                                <path
                                    d="M0.146447 0.146447C0.338683 -0.0478972 0.645911 -0.0270359 0.853553 0.146447L4.5 3.793L8.14645 0.146447C8.34171 -0.0488155 8.65829 -0.0488155 8.85355 0.146447C9.04882 0.341709 9.04882 0.658291 8.85355 0.853553L5.207 4.5L8.85355 8.14645C9.05934 8.35223 9.03129 8.67582 8.85355 8.85355C8.67582 9.03129 8.35409 9.02703 8.14645 8.85355L4.5 5.207L0.853553 8.85355C0.658291 9.04882 0.341709 9.04882 0.146447 8.85355C-0.0488155 8.65829 -0.0488155 8.34171 0.146447 8.14645L3.793 4.5L0.146447 0.853553C-0.0268697 0.680237 -0.0457894 0.34079 0.146447 0.146447Z"
                                />
                            </svg>
                        </slot>
                    </button>
                </slot>
            </div>
            ${n(t,e)}
        </div>
    </template>
`;class Ys extends rt{}class Xs extends(he(Ys)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class Qs extends Xs{readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.maxLength=this.maxlength,this.validate())}minlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.minLength=this.minlength,this.validate())}patternChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.pattern=this.pattern,this.validate())}sizeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.size=this.size)}spellcheckChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.spellcheck=this.spellcheck)}connectedCallback(){super.connectedCallback(),this.validate(),this.autofocus&&s.DOM.queueUpdate((()=>{this.focus()}))}validate(){super.validate(this.control)}handleTextInput(){this.value=this.control.value}handleClearInput(){this.value="",this.control.focus(),this.handleChange()}handleChange(){this.$emit("change")}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Qs.prototype,"readOnly",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Qs.prototype,"autofocus",void 0),(0,d.gn)([s.attr],Qs.prototype,"placeholder",void 0),(0,d.gn)([s.attr],Qs.prototype,"list",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Qs.prototype,"maxlength",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Qs.prototype,"minlength",void 0),(0,d.gn)([s.attr],Qs.prototype,"pattern",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],Qs.prototype,"size",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Qs.prototype,"spellcheck",void 0),(0,d.gn)([s.observable],Qs.prototype,"defaultSlottedNodes",void 0);class Zs{}dt(Zs,Rt),dt(Qs,o,Zs);class Js extends os{}class to extends(he(Js)){constructor(){super(...arguments),this.proxy=document.createElement("select")}}class eo extends to{constructor(){super(...arguments),this.open=!1,this.forcedPosition=!1,this.listboxId=Ne("listbox-"),this.maxHeight=0}openChanged(t,e){if(this.collapsible){if(this.open)return this.ariaControls=this.listboxId,this.ariaExpanded="true",this.setPositioning(),this.focusAndScrollOptionIntoView(),this.indexWhenOpened=this.selectedIndex,void s.DOM.queueUpdate((()=>this.focus()));this.ariaControls="",this.ariaExpanded="false"}}get collapsible(){return!(this.multiple||"number"==typeof this.size)}get value(){return s.Observable.track(this,"value"),this._value}set value(t){var e,i,o,n,a,r,l;const h=`${this._value}`;if(null===(e=this._options)||void 0===e?void 0:e.length){const e=this._options.findIndex((e=>e.value===t)),s=null!==(o=null===(i=this._options[this.selectedIndex])||void 0===i?void 0:i.value)&&void 0!==o?o:null,h=null!==(a=null===(n=this._options[e])||void 0===n?void 0:n.value)&&void 0!==a?a:null;-1!==e&&s===h||(t="",this.selectedIndex=e),t=null!==(l=null===(r=this.firstSelectedOption)||void 0===r?void 0:r.value)&&void 0!==l?l:t}h!==t&&(this._value=t,super.valueChanged(h,t),s.Observable.notify(this,"value"),this.updateDisplayValue())}updateValue(t){var e,i;this.$fastController.isConnected&&(this.value=null!==(i=null===(e=this.firstSelectedOption)||void 0===e?void 0:e.value)&&void 0!==i?i:""),t&&(this.$emit("input"),this.$emit("change",this,{bubbles:!0,composed:void 0}))}selectedIndexChanged(t,e){super.selectedIndexChanged(t,e),this.updateValue()}positionChanged(t,e){this.positionAttribute=e,this.setPositioning()}setPositioning(){const t=this.getBoundingClientRect(),e=window.innerHeight-t.bottom;this.position=this.forcedPosition?this.positionAttribute:t.top>e?We.above:We.below,this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position,this.maxHeight=this.position===We.above?~~t.top:~~e}get displayValue(){var t,e;return s.Observable.track(this,"displayValue"),null!==(e=null===(t=this.firstSelectedOption)||void 0===t?void 0:t.text)&&void 0!==e?e:""}disabledChanged(t,e){super.disabledChanged&&super.disabledChanged(t,e),this.ariaDisabled=this.disabled?"true":"false"}formResetCallback(){this.setProxyOptions(),super.setDefaultSelectedOption(),-1===this.selectedIndex&&(this.selectedIndex=0)}clickHandler(t){if(!this.disabled){if(this.open){const e=t.target.closest("option,[role=option]");if(e&&e.disabled)return}return super.clickHandler(t),this.open=this.collapsible&&!this.open,this.open||this.indexWhenOpened===this.selectedIndex||this.updateValue(!0),!0}}focusoutHandler(t){var e;if(super.focusoutHandler(t),!this.open)return!0;const i=t.relatedTarget;this.isSameNode(i)?this.focus():(null===(e=this.options)||void 0===e?void 0:e.includes(i))||(this.open=!1,this.indexWhenOpened!==this.selectedIndex&&this.updateValue(!0))}handleChange(t,e){super.handleChange(t,e),"value"===e&&this.updateValue()}slottedOptionsChanged(t,e){this.options.forEach((t=>{s.Observable.getNotifier(t).unsubscribe(this,"value")})),super.slottedOptionsChanged(t,e),this.options.forEach((t=>{s.Observable.getNotifier(t).subscribe(this,"value")})),this.setProxyOptions(),this.updateValue()}mousedownHandler(t){var e;return t.offsetX>=0&&t.offsetX<=(null===(e=this.listbox)||void 0===e?void 0:e.scrollWidth)?super.mousedownHandler(t):this.collapsible}multipleChanged(t,e){super.multipleChanged(t,e),this.proxy&&(this.proxy.multiple=e)}selectedOptionsChanged(t,e){var i;super.selectedOptionsChanged(t,e),null===(i=this.options)||void 0===i||i.forEach(((t,e)=>{var i;const s=null===(i=this.proxy)||void 0===i?void 0:i.options.item(e);s&&(s.selected=t.selected)}))}setDefaultSelectedOption(){var t;const e=null!==(t=this.options)&&void 0!==t?t:Array.from(this.children).filter(_e.slottedOptionFilter),i=null==e?void 0:e.findIndex((t=>t.hasAttribute("selected")||t.selected||t.value===this.value));this.selectedIndex=-1===i?0:i}setProxyOptions(){this.proxy instanceof HTMLSelectElement&&this.options&&(this.proxy.options.length=0,this.options.forEach((t=>{const e=t.proxy||(t instanceof HTMLOptionElement?t.cloneNode():null);e&&this.proxy.options.add(e)})))}keydownHandler(t){super.keydownHandler(t);const e=t.key||t.key.charCodeAt(0);switch(e){case wt:t.preventDefault(),this.collapsible&&this.typeAheadExpired&&(this.open=!this.open);break;case Ct:case xt:t.preventDefault();break;case ft:t.preventDefault(),this.open=!this.open;break;case yt:this.collapsible&&this.open&&(t.preventDefault(),this.open=!1);break;case $t:return this.collapsible&&this.open&&(t.preventDefault(),this.open=!1),!0}return this.open||this.indexWhenOpened===this.selectedIndex||(this.updateValue(!0),this.indexWhenOpened=this.selectedIndex),!(e===mt||e===bt)}connectedCallback(){super.connectedCallback(),this.forcedPosition=!!this.positionAttribute,this.addEventListener("contentchange",this.updateDisplayValue)}disconnectedCallback(){this.removeEventListener("contentchange",this.updateDisplayValue),super.disconnectedCallback()}sizeChanged(t,e){super.sizeChanged(t,e),this.proxy&&(this.proxy.size=e)}updateDisplayValue(){this.collapsible&&s.Observable.notify(this,"displayValue")}}(0,d.gn)([(0,s.attr)({attribute:"open",mode:"boolean"})],eo.prototype,"open",void 0),(0,d.gn)([s.volatile],eo.prototype,"collapsible",null),(0,d.gn)([s.observable],eo.prototype,"control",void 0),(0,d.gn)([(0,s.attr)({attribute:"position"})],eo.prototype,"positionAttribute",void 0),(0,d.gn)([s.observable],eo.prototype,"position",void 0),(0,d.gn)([s.observable],eo.prototype,"maxHeight",void 0);class io{}(0,d.gn)([s.observable],io.prototype,"ariaControls",void 0),dt(io,Ke),dt(eo,o,io);const so=(t,e)=>s.html`
    <template
        class="${t=>[t.collapsible&&"collapsible",t.collapsible&&t.open&&"open",t.disabled&&"disabled",t.collapsible&&t.position].filter(Boolean).join(" ")}"
        aria-activedescendant="${t=>t.ariaActiveDescendant}"
        aria-controls="${t=>t.ariaControls}"
        aria-disabled="${t=>t.ariaDisabled}"
        aria-expanded="${t=>t.ariaExpanded}"
        aria-haspopup="${t=>t.collapsible?"listbox":null}"
        aria-multiselectable="${t=>t.ariaMultiSelectable}"
        ?open="${t=>t.open}"
        role="combobox"
        tabindex="${t=>t.disabled?null:"0"}"
        @click="${(t,e)=>t.clickHandler(e.event)}"
        @focusin="${(t,e)=>t.focusinHandler(e.event)}"
        @focusout="${(t,e)=>t.focusoutHandler(e.event)}"
        @keydown="${(t,e)=>t.keydownHandler(e.event)}"
        @mousedown="${(t,e)=>t.mousedownHandler(e.event)}"
    >
        ${(0,s.when)((t=>t.collapsible),s.html`
                <div
                    class="control"
                    part="control"
                    ?disabled="${t=>t.disabled}"
                    ${(0,s.ref)("control")}
                >
                    ${a(t,e)}
                    <slot name="button-container">
                        <div class="selected-value" part="selected-value">
                            <slot name="selected-value">${t=>t.displayValue}</slot>
                        </div>
                        <div aria-hidden="true" class="indicator" part="indicator">
                            <slot name="indicator">
                                ${e.indicator||""}
                            </slot>
                        </div>
                    </slot>
                    ${n(t,e)}
                </div>
            `)}
        <div
            class="listbox"
            id="${t=>t.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${t=>t.disabled}"
            ?hidden="${t=>!!t.collapsible&&!t.open}"
            ${(0,s.ref)("listbox")}
        >
            <slot
                ${(0,s.slotted)({filter:_e.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`,oo=(t,e)=>s.html`
    <template
        class="${t=>"circle"===t.shape?"circle":"rect"}"
        pattern="${t=>t.pattern}"
        ?shimmer="${t=>t.shimmer}"
    >
        ${(0,s.when)((t=>!0===t.shimmer),s.html`
                <span class="shimmer"></span>
            `)}
        <object type="image/svg+xml" data="${t=>t.pattern}" role="presentation">
            <img class="pattern" src="${t=>t.pattern}" />
        </object>
        <slot></slot>
    </template>
`;class no extends rt{constructor(){super(...arguments),this.shape="rect"}}(0,d.gn)([s.attr],no.prototype,"fill",void 0),(0,d.gn)([s.attr],no.prototype,"shape",void 0),(0,d.gn)([s.attr],no.prototype,"pattern",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],no.prototype,"shimmer",void 0);const ao=(t,e)=>s.html`
    <template
        aria-disabled="${t=>t.disabled}"
        class="${t=>t.sliderOrientation||Qi.i.horizontal}
            ${t=>t.disabled?"disabled":""}"
    >
        <div ${(0,s.ref)("root")} part="root" class="root" style="${t=>t.positionStyle}">
            <div class="container">
                ${(0,s.when)((t=>!t.hideMark),s.html`
                        <div class="mark"></div>
                    `)}
                <div class="label">
                    <slot></slot>
                </div>
            </div>
        </div>
    </template>
`;function ro(t,e,i,s){let o=kt(0,1,(t-e)/(i-e));return s===Lt.N.rtl&&(o=1-o),o}const lo={min:0,max:0,direction:Lt.N.ltr,orientation:Qi.i.horizontal,disabled:!1};class ho extends rt{constructor(){super(...arguments),this.hideMark=!1,this.sliderDirection=Lt.N.ltr,this.getSliderConfiguration=()=>{if(this.isSliderConfig(this.parentNode)){const t=this.parentNode,{min:e,max:i,direction:s,orientation:o,disabled:n}=t;void 0!==n&&(this.disabled=n),this.sliderDirection=s||Lt.N.ltr,this.sliderOrientation=o||Qi.i.horizontal,this.sliderMaxPosition=i,this.sliderMinPosition=e}else this.sliderDirection=lo.direction||Lt.N.ltr,this.sliderOrientation=lo.orientation||Qi.i.horizontal,this.sliderMaxPosition=lo.max,this.sliderMinPosition=lo.min},this.positionAsStyle=()=>{const t=this.sliderDirection?this.sliderDirection:Lt.N.ltr,e=ro(Number(this.position),Number(this.sliderMinPosition),Number(this.sliderMaxPosition));let i=Math.round(100*(1-e)),s=Math.round(100*e);return Number.isNaN(s)&&Number.isNaN(i)&&(i=50,s=50),this.sliderOrientation===Qi.i.horizontal?t===Lt.N.rtl?`right: ${s}%; left: ${i}%;`:`left: ${s}%; right: ${i}%;`:`top: ${s}%; bottom: ${i}%;`}}positionChanged(){this.positionStyle=this.positionAsStyle()}sliderOrientationChanged(){}connectedCallback(){super.connectedCallback(),this.getSliderConfiguration(),this.positionStyle=this.positionAsStyle(),this.notifier=s.Observable.getNotifier(this.parentNode),this.notifier.subscribe(this,"orientation"),this.notifier.subscribe(this,"direction"),this.notifier.subscribe(this,"max"),this.notifier.subscribe(this,"min")}disconnectedCallback(){super.disconnectedCallback(),this.notifier.unsubscribe(this,"orientation"),this.notifier.unsubscribe(this,"direction"),this.notifier.unsubscribe(this,"max"),this.notifier.unsubscribe(this,"min")}handleChange(t,e){switch(e){case"direction":this.sliderDirection=t.direction;break;case"orientation":this.sliderOrientation=t.orientation;break;case"max":this.sliderMaxPosition=t.max;break;case"min":this.sliderMinPosition=t.min}this.positionStyle=this.positionAsStyle()}isSliderConfig(t){return void 0!==t.max&&void 0!==t.min}}(0,d.gn)([s.observable],ho.prototype,"positionStyle",void 0),(0,d.gn)([s.attr],ho.prototype,"position",void 0),(0,d.gn)([(0,s.attr)({attribute:"hide-mark",mode:"boolean"})],ho.prototype,"hideMark",void 0),(0,d.gn)([(0,s.attr)({attribute:"disabled",mode:"boolean"})],ho.prototype,"disabled",void 0),(0,d.gn)([s.observable],ho.prototype,"sliderOrientation",void 0),(0,d.gn)([s.observable],ho.prototype,"sliderMinPosition",void 0),(0,d.gn)([s.observable],ho.prototype,"sliderMaxPosition",void 0),(0,d.gn)([s.observable],ho.prototype,"sliderDirection",void 0);const co=(t,e)=>s.html`
    <template
        role="slider"
        class="${t=>t.readOnly?"readonly":""}
        ${t=>t.orientation||Qi.i.horizontal}"
        tabindex="${t=>t.disabled?null:0}"
        aria-valuetext="${t=>t.valueTextFormatter(t.value)}"
        aria-valuenow="${t=>t.value}"
        aria-valuemin="${t=>t.min}"
        aria-valuemax="${t=>t.max}"
        aria-disabled="${t=>!!t.disabled||void 0}"
        aria-readonly="${t=>!!t.readOnly||void 0}"
        aria-orientation="${t=>t.orientation}"
        class="${t=>t.orientation}"
    >
        <div part="positioning-region" class="positioning-region">
            <div ${(0,s.ref)("track")} part="track-container" class="track">
                <slot name="track"></slot>
                <div part="track-start" class="track-start" style="${t=>t.position}">
                    <slot name="track-start"></slot>
                </div>
            </div>
            <slot></slot>
            <div
                ${(0,s.ref)("thumb")}
                part="thumb-container"
                class="thumb-container"
                style="${t=>t.position}"
            >
                <slot name="thumb">${e.thumb||""}</slot>
            </div>
        </div>
    </template>
`;class uo extends rt{}class po extends(he(uo)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}const mo={singleValue:"single-value"};class vo extends po{constructor(){super(...arguments),this.direction=Lt.N.ltr,this.isDragging=!1,this.trackWidth=0,this.trackMinWidth=0,this.trackHeight=0,this.trackLeft=0,this.trackMinHeight=0,this.valueTextFormatter=()=>null,this.min=0,this.max=10,this.step=1,this.orientation=Qi.i.horizontal,this.mode=mo.singleValue,this.keypressHandler=t=>{if(!this.readOnly)if(t.key===Ct)t.preventDefault(),this.value=`${this.min}`;else if(t.key===xt)t.preventDefault(),this.value=`${this.max}`;else if(!t.shiftKey)switch(t.key){case gt:case bt:t.preventDefault(),this.increment();break;case vt:case mt:t.preventDefault(),this.decrement()}},this.setupTrackConstraints=()=>{const t=this.track.getBoundingClientRect();this.trackWidth=this.track.clientWidth,this.trackMinWidth=this.track.clientLeft,this.trackHeight=t.bottom,this.trackMinHeight=t.top,this.trackLeft=this.getBoundingClientRect().left,0===this.trackWidth&&(this.trackWidth=1)},this.setupListeners=(t=!1)=>{const e=(t?"remove":"add")+"EventListener";this[e]("keydown",this.keypressHandler),this[e]("mousedown",this.handleMouseDown),this.thumb[e]("mousedown",this.handleThumbMouseDown,{passive:!0}),this.thumb[e]("touchstart",this.handleThumbMouseDown,{passive:!0}),t&&(this.handleMouseDown(null),this.handleThumbMouseDown(null))},this.initialValue="",this.handleThumbMouseDown=t=>{if(t){if(this.readOnly||this.disabled||t.defaultPrevented)return;t.target.focus()}const e=(null!==t?"add":"remove")+"EventListener";window[e]("mouseup",this.handleWindowMouseUp),window[e]("mousemove",this.handleMouseMove,{passive:!0}),window[e]("touchmove",this.handleMouseMove,{passive:!0}),window[e]("touchend",this.handleWindowMouseUp),this.isDragging=null!==t},this.handleMouseMove=t=>{if(this.readOnly||this.disabled||t.defaultPrevented)return;const e=window.TouchEvent&&t instanceof TouchEvent?t.touches[0]:t,i=this.orientation===Qi.i.horizontal?e.pageX-document.documentElement.scrollLeft-this.trackLeft:e.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(i)}`},this.calculateNewValue=t=>{const e=ro(t,this.orientation===Qi.i.horizontal?this.trackMinWidth:this.trackMinHeight,this.orientation===Qi.i.horizontal?this.trackWidth:this.trackHeight,this.direction),i=(this.max-this.min)*e+this.min;return this.convertToConstrainedValue(i)},this.handleWindowMouseUp=t=>{this.stopDragging()},this.stopDragging=()=>{this.isDragging=!1,this.handleMouseDown(null),this.handleThumbMouseDown(null)},this.handleMouseDown=t=>{const e=(null!==t?"add":"remove")+"EventListener";if((null===t||!this.disabled&&!this.readOnly)&&(window[e]("mouseup",this.handleWindowMouseUp),window.document[e]("mouseleave",this.handleWindowMouseUp),window[e]("mousemove",this.handleMouseMove),t)){t.preventDefault(),this.setupTrackConstraints(),t.target.focus();const e=this.orientation===Qi.i.horizontal?t.pageX-document.documentElement.scrollLeft-this.trackLeft:t.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(e)}`}},this.convertToConstrainedValue=t=>{isNaN(t)&&(t=this.min);let e=t-this.min;const i=e-Math.round(e/this.step)*(this.stepMultiplier*this.step)/this.stepMultiplier;return e=i>=Number(this.step)/2?e-i+Number(this.step):e-i,e+this.min}}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(t){this.value=t.toString()}valueChanged(t,e){super.valueChanged(t,e),this.$fastController.isConnected&&this.setThumbPositionForOrientation(this.direction),this.$emit("change")}minChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.min=`${this.min}`),this.validate()}maxChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.max=`${this.max}`),this.validate()}stepChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.step=`${this.step}`),this.updateStepMultiplier(),this.validate()}orientationChanged(){this.$fastController.isConnected&&this.setThumbPositionForOrientation(this.direction)}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type","range"),this.direction=Bt(this),this.updateStepMultiplier(),this.setupTrackConstraints(),this.setupListeners(),this.setupDefaultValue(),this.setThumbPositionForOrientation(this.direction)}disconnectedCallback(){this.setupListeners(!0)}increment(){const t=this.direction!==Lt.N.rtl&&this.orientation!==Qi.i.vertical?Number(this.value)+Number(this.step):Number(this.value)-Number(this.step),e=this.convertToConstrainedValue(t),i=e<Number(this.max)?`${e}`:`${this.max}`;this.value=i}decrement(){const t=this.direction!==Lt.N.rtl&&this.orientation!==Qi.i.vertical?Number(this.value)-Number(this.step):Number(this.value)+Number(this.step),e=this.convertToConstrainedValue(t),i=e>Number(this.min)?`${e}`:`${this.min}`;this.value=i}setThumbPositionForOrientation(t){const e=100*(1-ro(Number(this.value),Number(this.min),Number(this.max),t));this.orientation===Qi.i.horizontal?this.position=this.isDragging?`right: ${e}%; transition: none;`:`right: ${e}%; transition: all 0.2s ease;`:this.position=this.isDragging?`bottom: ${e}%; transition: none;`:`bottom: ${e}%; transition: all 0.2s ease;`}updateStepMultiplier(){const t=this.step+"",e=this.step%1?t.length-t.indexOf(".")-1:0;this.stepMultiplier=Math.pow(10,e)}get midpoint(){return`${this.convertToConstrainedValue((this.max+this.min)/2)}`}setupDefaultValue(){if("string"==typeof this.value)if(0===this.value.length)this.initialValue=this.midpoint;else{const t=parseFloat(this.value);!Number.isNaN(t)&&(t<this.min||t>this.max)&&(this.value=this.midpoint)}}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],vo.prototype,"readOnly",void 0),(0,d.gn)([s.observable],vo.prototype,"direction",void 0),(0,d.gn)([s.observable],vo.prototype,"isDragging",void 0),(0,d.gn)([s.observable],vo.prototype,"position",void 0),(0,d.gn)([s.observable],vo.prototype,"trackWidth",void 0),(0,d.gn)([s.observable],vo.prototype,"trackMinWidth",void 0),(0,d.gn)([s.observable],vo.prototype,"trackHeight",void 0),(0,d.gn)([s.observable],vo.prototype,"trackLeft",void 0),(0,d.gn)([s.observable],vo.prototype,"trackMinHeight",void 0),(0,d.gn)([s.observable],vo.prototype,"valueTextFormatter",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],vo.prototype,"min",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],vo.prototype,"max",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],vo.prototype,"step",void 0),(0,d.gn)([s.attr],vo.prototype,"orientation",void 0),(0,d.gn)([s.attr],vo.prototype,"mode",void 0);const go=(t,e)=>s.html`
    <template
        role="switch"
        aria-checked="${t=>t.checked}"
        aria-disabled="${t=>t.disabled}"
        aria-readonly="${t=>t.readOnly}"
        tabindex="${t=>t.disabled?null:0}"
        @keypress="${(t,e)=>t.keypressHandler(e.event)}"
        @click="${(t,e)=>t.clickHandler(e.event)}"
        class="${t=>t.checked?"checked":""}"
    >
        <label
            part="label"
            class="${t=>t.defaultSlottedNodes&&t.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
        <div part="switch" class="switch">
            <slot name="switch">${e.switch||""}</slot>
        </div>
        <span class="status-message" part="status-message">
            <span class="checked-message" part="checked-message">
                <slot name="checked-message"></slot>
            </span>
            <span class="unchecked-message" part="unchecked-message">
                <slot name="unchecked-message"></slot>
            </span>
        </span>
    </template>
`;class bo extends rt{}class fo extends(de(bo)){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class yo extends fo{constructor(){super(),this.initialValue="on",this.keypressHandler=t=>{if(!this.readOnly)switch(t.key){case ft:case wt:this.checked=!this.checked}},this.clickHandler=t=>{this.disabled||this.readOnly||(this.checked=!this.checked)},this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly),this.readOnly?this.classList.add("readonly"):this.classList.remove("readonly")}checkedChanged(t,e){super.checkedChanged(t,e),this.checked?this.classList.add("checked"):this.classList.remove("checked")}}(0,d.gn)([(0,s.attr)({attribute:"readonly",mode:"boolean"})],yo.prototype,"readOnly",void 0),(0,d.gn)([s.observable],yo.prototype,"defaultSlottedNodes",void 0);const Co=(t,e)=>s.html`
    <template slot="tabpanel" role="tabpanel">
        <slot></slot>
    </template>
`;class xo extends rt{}const wo=(t,e)=>s.html`
    <template slot="tab" role="tab" aria-disabled="${t=>t.disabled}">
        <slot></slot>
    </template>
`;class $o extends rt{}(0,d.gn)([(0,s.attr)({mode:"boolean"})],$o.prototype,"disabled",void 0);const Io=(t,e)=>s.html`
    <template class="${t=>t.orientation}">
        ${a(t,e)}
        <div class="tablist" part="tablist" role="tablist">
            <slot class="tab" name="tab" part="tab" ${(0,s.slotted)("tabs")}></slot>

            ${(0,s.when)((t=>t.showActiveIndicator),s.html`
                    <div
                        ${(0,s.ref)("activeIndicatorRef")}
                        class="activeIndicator"
                        part="activeIndicator"
                    ></div>
                `)}
        </div>
        ${n(t,e)}
        <div class="tabpanel" part="tabpanel">
            <slot name="tabpanel" ${(0,s.slotted)("tabpanels")}></slot>
        </div>
    </template>
`,ko={vertical:"vertical",horizontal:"horizontal"};class Eo extends rt{constructor(){super(...arguments),this.orientation=ko.horizontal,this.activeindicator=!0,this.showActiveIndicator=!0,this.prevActiveTabIndex=0,this.activeTabIndex=0,this.ticking=!1,this.change=()=>{this.$emit("change",this.activetab)},this.isDisabledElement=t=>"true"===t.getAttribute("aria-disabled"),this.isHiddenElement=t=>t.hasAttribute("hidden"),this.isFocusableElement=t=>!this.isDisabledElement(t)&&!this.isHiddenElement(t),this.setTabs=()=>{const t="gridColumn",e="gridRow",i=this.isHorizontal()?t:e;this.activeTabIndex=this.getActiveIndex(),this.showActiveIndicator=!1,this.tabs.forEach(((s,o)=>{if("tab"===s.slot){const t=this.activeTabIndex===o&&this.isFocusableElement(s);this.activeindicator&&this.isFocusableElement(s)&&(this.showActiveIndicator=!0);const e=this.tabIds[o],i=this.tabpanelIds[o];s.setAttribute("id",e),s.setAttribute("aria-selected",t?"true":"false"),s.setAttribute("aria-controls",i),s.addEventListener("click",this.handleTabClick),s.addEventListener("keydown",this.handleTabKeyDown),s.setAttribute("tabindex",t?"0":"-1"),t&&(this.activetab=s,this.activeid=e)}s.style[t]="",s.style[e]="",s.style[i]=`${o+1}`,this.isHorizontal()?s.classList.remove("vertical"):s.classList.add("vertical")}))},this.setTabPanels=()=>{this.tabpanels.forEach(((t,e)=>{const i=this.tabIds[e],s=this.tabpanelIds[e];t.setAttribute("id",s),t.setAttribute("aria-labelledby",i),this.activeTabIndex!==e?t.setAttribute("hidden",""):t.removeAttribute("hidden")}))},this.handleTabClick=t=>{const e=t.currentTarget;1===e.nodeType&&this.isFocusableElement(e)&&(this.prevActiveTabIndex=this.activeTabIndex,this.activeTabIndex=this.tabs.indexOf(e),this.setComponent())},this.handleTabKeyDown=t=>{if(this.isHorizontal())switch(t.key){case vt:t.preventDefault(),this.adjustBackward(t);break;case gt:t.preventDefault(),this.adjustForward(t)}else switch(t.key){case bt:t.preventDefault(),this.adjustBackward(t);break;case mt:t.preventDefault(),this.adjustForward(t)}switch(t.key){case Ct:t.preventDefault(),this.adjust(-this.activeTabIndex);break;case xt:t.preventDefault(),this.adjust(this.tabs.length-this.activeTabIndex-1)}},this.adjustForward=t=>{const e=this.tabs;let i=0;for(i=this.activetab?e.indexOf(this.activetab)+1:1,i===e.length&&(i=0);i<e.length&&e.length>1;){if(this.isFocusableElement(e[i])){this.moveToTabByIndex(e,i);break}if(this.activetab&&i===e.indexOf(this.activetab))break;i+1>=e.length?i=0:i+=1}},this.adjustBackward=t=>{const e=this.tabs;let i=0;for(i=this.activetab?e.indexOf(this.activetab)-1:0,i=i<0?e.length-1:i;i>=0&&e.length>1;){if(this.isFocusableElement(e[i])){this.moveToTabByIndex(e,i);break}i-1<0?i=e.length-1:i-=1}},this.moveToTabByIndex=(t,e)=>{const i=t[e];this.activetab=i,this.prevActiveTabIndex=this.activeTabIndex,this.activeTabIndex=e,i.focus(),this.setComponent()}}orientationChanged(){this.$fastController.isConnected&&(this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}activeidChanged(t,e){this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length&&(this.prevActiveTabIndex=this.tabs.findIndex((e=>e.id===t)),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}tabsChanged(){this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length&&(this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}tabpanelsChanged(){this.$fastController.isConnected&&this.tabpanels.length<=this.tabs.length&&(this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}getActiveIndex(){return void 0!==this.activeid?-1===this.tabIds.indexOf(this.activeid)?0:this.tabIds.indexOf(this.activeid):0}getTabIds(){return this.tabs.map((t=>{var e;return null!==(e=t.getAttribute("id"))&&void 0!==e?e:`tab-${Ne()}`}))}getTabPanelIds(){return this.tabpanels.map((t=>{var e;return null!==(e=t.getAttribute("id"))&&void 0!==e?e:`panel-${Ne()}`}))}setComponent(){this.activeTabIndex!==this.prevActiveTabIndex&&(this.activeid=this.tabIds[this.activeTabIndex],this.focusTab(),this.change())}isHorizontal(){return this.orientation===ko.horizontal}handleActiveIndicatorPosition(){this.showActiveIndicator&&this.activeindicator&&this.activeTabIndex!==this.prevActiveTabIndex&&(this.ticking?this.ticking=!1:(this.ticking=!0,this.animateActiveIndicator()))}animateActiveIndicator(){this.ticking=!0;const t=this.isHorizontal()?"gridColumn":"gridRow",e=this.isHorizontal()?"translateX":"translateY",i=this.isHorizontal()?"offsetLeft":"offsetTop",s=this.activeIndicatorRef[i];this.activeIndicatorRef.style[t]=`${this.activeTabIndex+1}`;const o=this.activeIndicatorRef[i];this.activeIndicatorRef.style[t]=`${this.prevActiveTabIndex+1}`;const n=o-s;this.activeIndicatorRef.style.transform=`${e}(${n}px)`,this.activeIndicatorRef.classList.add("activeIndicatorTransition"),this.activeIndicatorRef.addEventListener("transitionend",(()=>{this.ticking=!1,this.activeIndicatorRef.style[t]=`${this.activeTabIndex+1}`,this.activeIndicatorRef.style.transform=`${e}(0px)`,this.activeIndicatorRef.classList.remove("activeIndicatorTransition")}))}adjust(t){const e=this.tabs.filter((t=>this.isFocusableElement(t))),i=e.indexOf(this.activetab),s=kt(0,e.length-1,i+t),o=this.tabs.indexOf(e[s]);o>-1&&this.moveToTabByIndex(this.tabs,o)}focusTab(){this.tabs[this.activeTabIndex].focus()}connectedCallback(){super.connectedCallback(),this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.activeTabIndex=this.getActiveIndex()}}(0,d.gn)([s.attr],Eo.prototype,"orientation",void 0),(0,d.gn)([s.attr],Eo.prototype,"activeid",void 0),(0,d.gn)([s.observable],Eo.prototype,"tabs",void 0),(0,d.gn)([s.observable],Eo.prototype,"tabpanels",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],Eo.prototype,"activeindicator",void 0),(0,d.gn)([s.observable],Eo.prototype,"activeIndicatorRef",void 0),(0,d.gn)([s.observable],Eo.prototype,"showActiveIndicator",void 0),dt(Eo,o);const To={none:"none",both:"both",horizontal:"horizontal",vertical:"vertical"},Oo=(t,e)=>s.html`
    <template
        class="
            ${t=>t.readOnly?"readonly":""}
            ${t=>t.resize!==To.none?`resize-${t.resize}`:""}"
    >
        <label
            part="label"
            for="control"
            class="${t=>t.defaultSlottedNodes&&t.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
        <textarea
            part="control"
            class="control"
            id="control"
            ?autofocus="${t=>t.autofocus}"
            cols="${t=>t.cols}"
            ?disabled="${t=>t.disabled}"
            form="${t=>t.form}"
            list="${t=>t.list}"
            maxlength="${t=>t.maxlength}"
            minlength="${t=>t.minlength}"
            name="${t=>t.name}"
            placeholder="${t=>t.placeholder}"
            ?readonly="${t=>t.readOnly}"
            ?required="${t=>t.required}"
            rows="${t=>t.rows}"
            ?spellcheck="${t=>t.spellcheck}"
            :value="${t=>t.value}"
            aria-atomic="${t=>t.ariaAtomic}"
            aria-busy="${t=>t.ariaBusy}"
            aria-controls="${t=>t.ariaControls}"
            aria-current="${t=>t.ariaCurrent}"
            aria-describedby="${t=>t.ariaDescribedby}"
            aria-details="${t=>t.ariaDetails}"
            aria-disabled="${t=>t.ariaDisabled}"
            aria-errormessage="${t=>t.ariaErrormessage}"
            aria-flowto="${t=>t.ariaFlowto}"
            aria-haspopup="${t=>t.ariaHaspopup}"
            aria-hidden="${t=>t.ariaHidden}"
            aria-invalid="${t=>t.ariaInvalid}"
            aria-keyshortcuts="${t=>t.ariaKeyshortcuts}"
            aria-label="${t=>t.ariaLabel}"
            aria-labelledby="${t=>t.ariaLabelledby}"
            aria-live="${t=>t.ariaLive}"
            aria-owns="${t=>t.ariaOwns}"
            aria-relevant="${t=>t.ariaRelevant}"
            aria-roledescription="${t=>t.ariaRoledescription}"
            @input="${(t,e)=>t.handleTextInput()}"
            @change="${t=>t.handleChange()}"
            ${(0,s.ref)("control")}
        ></textarea>
    </template>
`;class Do extends rt{}class Ro extends(he(Do)){constructor(){super(...arguments),this.proxy=document.createElement("textarea")}}class So extends Ro{constructor(){super(...arguments),this.resize=To.none,this.cols=20,this.handleTextInput=()=>{this.value=this.control.value}}readOnlyChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.readOnly=this.readOnly)}autofocusChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.autofocus=this.autofocus)}listChanged(){this.proxy instanceof HTMLTextAreaElement&&this.proxy.setAttribute("list",this.list)}maxlengthChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.maxLength=this.maxlength)}minlengthChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.minLength=this.minlength)}spellcheckChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.spellcheck=this.spellcheck)}select(){this.control.select(),this.$emit("select")}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],So.prototype,"readOnly",void 0),(0,d.gn)([s.attr],So.prototype,"resize",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],So.prototype,"autofocus",void 0),(0,d.gn)([(0,s.attr)({attribute:"form"})],So.prototype,"formId",void 0),(0,d.gn)([s.attr],So.prototype,"list",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],So.prototype,"maxlength",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter})],So.prototype,"minlength",void 0),(0,d.gn)([s.attr],So.prototype,"name",void 0),(0,d.gn)([s.attr],So.prototype,"placeholder",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter,mode:"fromView"})],So.prototype,"cols",void 0),(0,d.gn)([(0,s.attr)({converter:s.nullableNumberConverter,mode:"fromView"})],So.prototype,"rows",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],So.prototype,"spellcheck",void 0),(0,d.gn)([s.observable],So.prototype,"defaultSlottedNodes",void 0),dt(So,As);const Ao=(t,e)=>s.html`
    <template
        class="
            ${t=>t.readOnly?"readonly":""}
        "
    >
        <label
            part="label"
            for="control"
            class="${t=>t.defaultSlottedNodes&&t.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot
                ${(0,s.slotted)({property:"defaultSlottedNodes",filter:Ws})}
            ></slot>
        </label>
        <div class="root" part="root">
            ${a(t,e)}
            <input
                class="control"
                part="control"
                id="control"
                @input="${t=>t.handleTextInput()}"
                @change="${t=>t.handleChange()}"
                ?autofocus="${t=>t.autofocus}"
                ?disabled="${t=>t.disabled}"
                list="${t=>t.list}"
                maxlength="${t=>t.maxlength}"
                minlength="${t=>t.minlength}"
                pattern="${t=>t.pattern}"
                placeholder="${t=>t.placeholder}"
                ?readonly="${t=>t.readOnly}"
                ?required="${t=>t.required}"
                size="${t=>t.size}"
                ?spellcheck="${t=>t.spellcheck}"
                :value="${t=>t.value}"
                type="${t=>t.type}"
                aria-atomic="${t=>t.ariaAtomic}"
                aria-busy="${t=>t.ariaBusy}"
                aria-controls="${t=>t.ariaControls}"
                aria-current="${t=>t.ariaCurrent}"
                aria-describedby="${t=>t.ariaDescribedby}"
                aria-details="${t=>t.ariaDetails}"
                aria-disabled="${t=>t.ariaDisabled}"
                aria-errormessage="${t=>t.ariaErrormessage}"
                aria-flowto="${t=>t.ariaFlowto}"
                aria-haspopup="${t=>t.ariaHaspopup}"
                aria-hidden="${t=>t.ariaHidden}"
                aria-invalid="${t=>t.ariaInvalid}"
                aria-keyshortcuts="${t=>t.ariaKeyshortcuts}"
                aria-label="${t=>t.ariaLabel}"
                aria-labelledby="${t=>t.ariaLabelledby}"
                aria-live="${t=>t.ariaLive}"
                aria-owns="${t=>t.ariaOwns}"
                aria-relevant="${t=>t.ariaRelevant}"
                aria-roledescription="${t=>t.ariaRoledescription}"
                ${(0,s.ref)("control")}
            />
            ${n(t,e)}
        </div>
    </template>
`,Fo=(t,e)=>s.html`
    <template
        aria-label="${t=>t.ariaLabel}"
        aria-labelledby="${t=>t.ariaLabelledby}"
        aria-orientation="${t=>t.orientation}"
        orientation="${t=>t.orientation}"
        role="toolbar"
        @mousedown="${(t,e)=>t.mouseDownHandler(e.event)}"
        @focusin="${(t,e)=>t.focusinHandler(e.event)}"
        @keydown="${(t,e)=>t.keydownHandler(e.event)}"
        ${(0,s.children)({property:"childItems",attributeFilter:["disabled","hidden"],filter:(0,s.elements)(),subtree:!0})}
    >
        <slot name="label"></slot>
        <div class="positioning-region" part="positioning-region">
            ${a(t,e)}
            <slot
                ${(0,s.slotted)({filter:(0,s.elements)(),property:"slottedItems"})}
            ></slot>
            ${n(t,e)}
        </div>
    </template>
`,Lo=Object.freeze({[It.ArrowUp]:{[Qi.i.vertical]:-1},[It.ArrowDown]:{[Qi.i.vertical]:1},[It.ArrowLeft]:{[Qi.i.horizontal]:{[Lt.N.ltr]:-1,[Lt.N.rtl]:1}},[It.ArrowRight]:{[Qi.i.horizontal]:{[Lt.N.ltr]:1,[Lt.N.rtl]:-1}}});class Po extends rt{constructor(){super(...arguments),this._activeIndex=0,this.direction=Lt.N.ltr,this.orientation=Qi.i.horizontal}get activeIndex(){return s.Observable.track(this,"activeIndex"),this._activeIndex}set activeIndex(t){this.$fastController.isConnected&&(this._activeIndex=kt(0,this.focusableElements.length-1,t),s.Observable.notify(this,"activeIndex"))}slottedItemsChanged(){this.$fastController.isConnected&&this.reduceFocusableElements()}mouseDownHandler(t){var e;const i=null===(e=this.focusableElements)||void 0===e?void 0:e.findIndex((e=>e.contains(t.target)));return i>-1&&this.activeIndex!==i&&this.setFocusedElement(i),!0}childItemsChanged(t,e){this.$fastController.isConnected&&this.reduceFocusableElements()}connectedCallback(){super.connectedCallback(),this.direction=Bt(this)}focusinHandler(t){const e=t.relatedTarget;e&&!this.contains(e)&&this.setFocusedElement()}getDirectionalIncrementer(t){var e,i,s,o,n;return null!==(n=null!==(s=null===(i=null===(e=Lo[t])||void 0===e?void 0:e[this.orientation])||void 0===i?void 0:i[this.direction])&&void 0!==s?s:null===(o=Lo[t])||void 0===o?void 0:o[this.orientation])&&void 0!==n?n:0}keydownHandler(t){const e=t.key;if(!(e in It)||t.defaultPrevented||t.shiftKey)return!0;const i=this.getDirectionalIncrementer(e);if(!i)return!t.target.closest("[role=radiogroup]");const s=this.activeIndex+i;return this.focusableElements[s]&&t.preventDefault(),this.setFocusedElement(s),!0}get allSlottedItems(){return[...this.start.assignedElements(),...this.slottedItems,...this.end.assignedElements()]}reduceFocusableElements(){var t;const e=null===(t=this.focusableElements)||void 0===t?void 0:t[this.activeIndex];this.focusableElements=this.allSlottedItems.reduce(Po.reduceFocusableItems,[]);const i=this.focusableElements.indexOf(e);this.activeIndex=Math.max(0,i),this.setFocusableElements()}setFocusedElement(t=this.activeIndex){var e;this.activeIndex=t,this.setFocusableElements(),null===(e=this.focusableElements[this.activeIndex])||void 0===e||e.focus()}static reduceFocusableItems(t,e){var i,s,o,n;const a="radio"===e.getAttribute("role"),r=null===(s=null===(i=e.$fastController)||void 0===i?void 0:i.definition.shadowOptions)||void 0===s?void 0:s.delegatesFocus,l=Array.from(null!==(n=null===(o=e.shadowRoot)||void 0===o?void 0:o.querySelectorAll("*"))&&void 0!==n?n:[]).some((t=>qi(t)));return e.hasAttribute("disabled")||e.hasAttribute("hidden")||!(qi(e)||a||r||l)?e.childElementCount?t.concat(Array.from(e.children).reduce(Po.reduceFocusableItems,[])):t:(t.push(e),t)}setFocusableElements(){this.$fastController.isConnected&&this.focusableElements.length>0&&this.focusableElements.forEach(((t,e)=>{t.tabIndex=this.activeIndex===e?0:-1}))}}(0,d.gn)([s.observable],Po.prototype,"direction",void 0),(0,d.gn)([s.attr],Po.prototype,"orientation",void 0),(0,d.gn)([s.observable],Po.prototype,"slottedItems",void 0),(0,d.gn)([s.observable],Po.prototype,"slottedLabel",void 0),(0,d.gn)([s.observable],Po.prototype,"childItems",void 0);class Mo{}(0,d.gn)([(0,s.attr)({attribute:"aria-labelledby"})],Mo.prototype,"ariaLabelledby",void 0),(0,d.gn)([(0,s.attr)({attribute:"aria-label"})],Mo.prototype,"ariaLabel",void 0),dt(Mo,Rt),dt(Po,o,Mo);const Ho=(t,e)=>s.html`
        ${(0,s.when)((t=>t.tooltipVisible),s.html`
            <${t.tagFor(qt)}
                fixed-placement="true"
                auto-update-mode="${t=>t.autoUpdateMode}"
                vertical-positioning-mode="${t=>t.verticalPositioningMode}"
                vertical-default-position="${t=>t.verticalDefaultPosition}"
                vertical-inset="${t=>t.verticalInset}"
                vertical-scaling="${t=>t.verticalScaling}"
                horizontal-positioning-mode="${t=>t.horizontalPositioningMode}"
                horizontal-default-position="${t=>t.horizontalDefaultPosition}"
                horizontal-scaling="${t=>t.horizontalScaling}"
                horizontal-inset="${t=>t.horizontalInset}"
                vertical-viewport-lock="${t=>t.horizontalViewportLock}"
                horizontal-viewport-lock="${t=>t.verticalViewportLock}"
                dir="${t=>t.currentDirection}"
                ${(0,s.ref)("region")}
            >
                <div class="tooltip" part="tooltip" role="tooltip">
                    <slot></slot>
                </div>
            </${t.tagFor(qt)}>
        `)}
    `,zo={top:"top",right:"right",bottom:"bottom",left:"left",start:"start",end:"end",topLeft:"top-left",topRight:"top-right",bottomLeft:"bottom-left",bottomRight:"bottom-right",topStart:"top-start",topEnd:"top-end",bottomStart:"bottom-start",bottomEnd:"bottom-end"};class Vo extends rt{constructor(){super(...arguments),this.anchor="",this.delay=300,this.autoUpdateMode="anchor",this.anchorElement=null,this.viewportElement=null,this.verticalPositioningMode="dynamic",this.horizontalPositioningMode="dynamic",this.horizontalInset="false",this.verticalInset="false",this.horizontalScaling="content",this.verticalScaling="content",this.verticalDefaultPosition=void 0,this.horizontalDefaultPosition=void 0,this.tooltipVisible=!1,this.currentDirection=Lt.N.ltr,this.showDelayTimer=null,this.hideDelayTimer=null,this.isAnchorHoveredFocused=!1,this.isRegionHovered=!1,this.handlePositionChange=t=>{this.classList.toggle("top","start"===this.region.verticalPosition),this.classList.toggle("bottom","end"===this.region.verticalPosition),this.classList.toggle("inset-top","insetStart"===this.region.verticalPosition),this.classList.toggle("inset-bottom","insetEnd"===this.region.verticalPosition),this.classList.toggle("center-vertical","center"===this.region.verticalPosition),this.classList.toggle("left","start"===this.region.horizontalPosition),this.classList.toggle("right","end"===this.region.horizontalPosition),this.classList.toggle("inset-left","insetStart"===this.region.horizontalPosition),this.classList.toggle("inset-right","insetEnd"===this.region.horizontalPosition),this.classList.toggle("center-horizontal","center"===this.region.horizontalPosition)},this.handleRegionMouseOver=t=>{this.isRegionHovered=!0},this.handleRegionMouseOut=t=>{this.isRegionHovered=!1,this.startHideDelayTimer()},this.handleAnchorMouseOver=t=>{this.tooltipVisible?this.isAnchorHoveredFocused=!0:this.startShowDelayTimer()},this.handleAnchorMouseOut=t=>{this.isAnchorHoveredFocused=!1,this.clearShowDelayTimer(),this.startHideDelayTimer()},this.handleAnchorFocusIn=t=>{this.startShowDelayTimer()},this.handleAnchorFocusOut=t=>{this.isAnchorHoveredFocused=!1,this.clearShowDelayTimer(),this.startHideDelayTimer()},this.startHideDelayTimer=()=>{this.clearHideDelayTimer(),this.tooltipVisible&&(this.hideDelayTimer=window.setTimeout((()=>{this.updateTooltipVisibility()}),60))},this.clearHideDelayTimer=()=>{null!==this.hideDelayTimer&&(clearTimeout(this.hideDelayTimer),this.hideDelayTimer=null)},this.startShowDelayTimer=()=>{this.isAnchorHoveredFocused||(this.delay>1?null===this.showDelayTimer&&(this.showDelayTimer=window.setTimeout((()=>{this.startHover()}),this.delay)):this.startHover())},this.startHover=()=>{this.isAnchorHoveredFocused=!0,this.updateTooltipVisibility()},this.clearShowDelayTimer=()=>{null!==this.showDelayTimer&&(clearTimeout(this.showDelayTimer),this.showDelayTimer=null)},this.getAnchor=()=>{const t=this.getRootNode();return t instanceof ShadowRoot?t.getElementById(this.anchor):document.getElementById(this.anchor)},this.handleDocumentKeydown=t=>{!t.defaultPrevented&&this.tooltipVisible&&t.key===yt&&(this.isAnchorHoveredFocused=!1,this.updateTooltipVisibility(),this.$emit("dismiss"))},this.updateTooltipVisibility=()=>{if(!1===this.visible)this.hideTooltip();else{if(!0===this.visible)return void this.showTooltip();if(this.isAnchorHoveredFocused||this.isRegionHovered)return void this.showTooltip();this.hideTooltip()}},this.showTooltip=()=>{this.tooltipVisible||(this.currentDirection=Bt(this),this.tooltipVisible=!0,document.addEventListener("keydown",this.handleDocumentKeydown),s.DOM.queueUpdate(this.setRegionProps))},this.hideTooltip=()=>{this.tooltipVisible&&(this.clearHideDelayTimer(),null!==this.region&&void 0!==this.region&&(this.region.removeEventListener("positionchange",this.handlePositionChange),this.region.viewportElement=null,this.region.anchorElement=null,this.region.removeEventListener("mouseover",this.handleRegionMouseOver),this.region.removeEventListener("mouseout",this.handleRegionMouseOut)),document.removeEventListener("keydown",this.handleDocumentKeydown),this.tooltipVisible=!1)},this.setRegionProps=()=>{this.tooltipVisible&&(this.region.viewportElement=this.viewportElement,this.region.anchorElement=this.anchorElement,this.region.addEventListener("positionchange",this.handlePositionChange),this.region.addEventListener("mouseover",this.handleRegionMouseOver,{passive:!0}),this.region.addEventListener("mouseout",this.handleRegionMouseOut,{passive:!0}))}}visibleChanged(){this.$fastController.isConnected&&(this.updateTooltipVisibility(),this.updateLayout())}anchorChanged(){this.$fastController.isConnected&&(this.anchorElement=this.getAnchor())}positionChanged(){this.$fastController.isConnected&&this.updateLayout()}anchorElementChanged(t){if(this.$fastController.isConnected){if(null!=t&&(t.removeEventListener("mouseover",this.handleAnchorMouseOver),t.removeEventListener("mouseout",this.handleAnchorMouseOut),t.removeEventListener("focusin",this.handleAnchorFocusIn),t.removeEventListener("focusout",this.handleAnchorFocusOut)),null!==this.anchorElement&&void 0!==this.anchorElement){this.anchorElement.addEventListener("mouseover",this.handleAnchorMouseOver,{passive:!0}),this.anchorElement.addEventListener("mouseout",this.handleAnchorMouseOut,{passive:!0}),this.anchorElement.addEventListener("focusin",this.handleAnchorFocusIn,{passive:!0}),this.anchorElement.addEventListener("focusout",this.handleAnchorFocusOut,{passive:!0});const t=this.anchorElement.id;null!==this.anchorElement.parentElement&&this.anchorElement.parentElement.querySelectorAll(":hover").forEach((e=>{e.id===t&&this.startShowDelayTimer()}))}null!==this.region&&void 0!==this.region&&this.tooltipVisible&&(this.region.anchorElement=this.anchorElement),this.updateLayout()}}viewportElementChanged(){null!==this.region&&void 0!==this.region&&(this.region.viewportElement=this.viewportElement),this.updateLayout()}connectedCallback(){super.connectedCallback(),this.anchorElement=this.getAnchor(),this.updateTooltipVisibility()}disconnectedCallback(){this.hideTooltip(),this.clearShowDelayTimer(),this.clearHideDelayTimer(),super.disconnectedCallback()}updateLayout(){switch(this.verticalPositioningMode="locktodefault",this.horizontalPositioningMode="locktodefault",this.position){case zo.top:case zo.bottom:this.verticalDefaultPosition=this.position,this.horizontalDefaultPosition="center";break;case zo.right:case zo.left:case zo.start:case zo.end:this.verticalDefaultPosition="center",this.horizontalDefaultPosition=this.position;break;case zo.topLeft:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="left";break;case zo.topRight:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="right";break;case zo.bottomLeft:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="left";break;case zo.bottomRight:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="right";break;case zo.topStart:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="start";break;case zo.topEnd:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="end";break;case zo.bottomStart:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="start";break;case zo.bottomEnd:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="end";break;default:this.verticalPositioningMode="dynamic",this.horizontalPositioningMode="dynamic",this.verticalDefaultPosition=void 0,this.horizontalDefaultPosition="center"}}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],Vo.prototype,"visible",void 0),(0,d.gn)([s.attr],Vo.prototype,"anchor",void 0),(0,d.gn)([s.attr],Vo.prototype,"delay",void 0),(0,d.gn)([s.attr],Vo.prototype,"position",void 0),(0,d.gn)([(0,s.attr)({attribute:"auto-update-mode"})],Vo.prototype,"autoUpdateMode",void 0),(0,d.gn)([(0,s.attr)({attribute:"horizontal-viewport-lock"})],Vo.prototype,"horizontalViewportLock",void 0),(0,d.gn)([(0,s.attr)({attribute:"vertical-viewport-lock"})],Vo.prototype,"verticalViewportLock",void 0),(0,d.gn)([s.observable],Vo.prototype,"anchorElement",void 0),(0,d.gn)([s.observable],Vo.prototype,"viewportElement",void 0),(0,d.gn)([s.observable],Vo.prototype,"verticalPositioningMode",void 0),(0,d.gn)([s.observable],Vo.prototype,"horizontalPositioningMode",void 0),(0,d.gn)([s.observable],Vo.prototype,"horizontalInset",void 0),(0,d.gn)([s.observable],Vo.prototype,"verticalInset",void 0),(0,d.gn)([s.observable],Vo.prototype,"horizontalScaling",void 0),(0,d.gn)([s.observable],Vo.prototype,"verticalScaling",void 0),(0,d.gn)([s.observable],Vo.prototype,"verticalDefaultPosition",void 0),(0,d.gn)([s.observable],Vo.prototype,"horizontalDefaultPosition",void 0),(0,d.gn)([s.observable],Vo.prototype,"tooltipVisible",void 0),(0,d.gn)([s.observable],Vo.prototype,"currentDirection",void 0);const No=(t,e)=>s.html`
    <template
        role="treeitem"
        slot="${t=>t.isNestedItem()?"item":void 0}"
        tabindex="-1"
        class="${t=>t.expanded?"expanded":""} ${t=>t.selected?"selected":""} ${t=>t.nested?"nested":""}
            ${t=>t.disabled?"disabled":""}"
        aria-expanded="${t=>t.childItems&&t.childItemLength()>0?t.expanded:void 0}"
        aria-selected="${t=>t.selected}"
        aria-disabled="${t=>t.disabled}"
        @focusin="${(t,e)=>t.handleFocus(e.event)}"
        @focusout="${(t,e)=>t.handleBlur(e.event)}"
        ${(0,s.children)({property:"childItems",filter:(0,s.elements)()})}
    >
        <div class="positioning-region" part="positioning-region">
            <div class="content-region" part="content-region">
                ${(0,s.when)((t=>t.childItems&&t.childItemLength()>0),s.html`
                        <div
                            aria-hidden="true"
                            class="expand-collapse-button"
                            part="expand-collapse-button"
                            @click="${(t,e)=>t.handleExpandCollapseButtonClick(e.event)}"
                            ${(0,s.ref)("expandCollapseButton")}
                        >
                            <slot name="expand-collapse-glyph">
                                ${e.expandCollapseGlyph||""}
                            </slot>
                        </div>
                    `)}
                ${a(t,e)}
                <slot></slot>
                ${n(t,e)}
            </div>
        </div>
        ${(0,s.when)((t=>t.childItems&&t.childItemLength()>0&&(t.expanded||t.renderCollapsedChildren)),s.html`
                <div role="group" class="items" part="items">
                    <slot name="item" ${(0,s.slotted)("items")}></slot>
                </div>
            `)}
    </template>
`;function Bo(t){return Be(t)&&"treeitem"===t.getAttribute("role")}class qo extends rt{constructor(){super(...arguments),this.expanded=!1,this.focusable=!1,this.isNestedItem=()=>Bo(this.parentElement),this.handleExpandCollapseButtonClick=t=>{this.disabled||t.defaultPrevented||(this.expanded=!this.expanded)},this.handleFocus=t=>{this.setAttribute("tabindex","0")},this.handleBlur=t=>{this.setAttribute("tabindex","-1")}}expandedChanged(){this.$fastController.isConnected&&this.$emit("expanded-change",this)}selectedChanged(){this.$fastController.isConnected&&this.$emit("selected-change",this)}itemsChanged(t,e){this.$fastController.isConnected&&this.items.forEach((t=>{Bo(t)&&(t.nested=!0)}))}static focusItem(t){t.focusable=!0,t.focus()}childItemLength(){const t=this.childItems.filter((t=>Bo(t)));return t?t.length:0}}(0,d.gn)([(0,s.attr)({mode:"boolean"})],qo.prototype,"expanded",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],qo.prototype,"selected",void 0),(0,d.gn)([(0,s.attr)({mode:"boolean"})],qo.prototype,"disabled",void 0),(0,d.gn)([s.observable],qo.prototype,"focusable",void 0),(0,d.gn)([s.observable],qo.prototype,"childItems",void 0),(0,d.gn)([s.observable],qo.prototype,"items",void 0),(0,d.gn)([s.observable],qo.prototype,"nested",void 0),(0,d.gn)([s.observable],qo.prototype,"renderCollapsedChildren",void 0),dt(qo,o);const Uo=(t,e)=>s.html`
    <template
        role="tree"
        ${(0,s.ref)("treeView")}
        @keydown="${(t,e)=>t.handleKeyDown(e.event)}"
        @focusin="${(t,e)=>t.handleFocus(e.event)}"
        @focusout="${(t,e)=>t.handleBlur(e.event)}"
        @click="${(t,e)=>t.handleClick(e.event)}"
        @selected-change="${(t,e)=>t.handleSelectedChange(e.event)}"
    >
        <slot ${(0,s.slotted)("slottedTreeItems")}></slot>
    </template>
`;class jo extends rt{constructor(){super(...arguments),this.currentFocused=null,this.handleFocus=t=>{if(!(this.slottedTreeItems.length<1))return t.target===this?(null===this.currentFocused&&(this.currentFocused=this.getValidFocusableItem()),void(null!==this.currentFocused&&qo.focusItem(this.currentFocused))):void(this.contains(t.target)&&(this.setAttribute("tabindex","-1"),this.currentFocused=t.target))},this.handleBlur=t=>{t.target instanceof HTMLElement&&(null===t.relatedTarget||!this.contains(t.relatedTarget))&&this.setAttribute("tabindex","0")},this.handleKeyDown=t=>{if(t.defaultPrevented)return;if(this.slottedTreeItems.length<1)return!0;const e=this.getVisibleNodes();switch(t.key){case Ct:return void(e.length&&qo.focusItem(e[0]));case xt:return void(e.length&&qo.focusItem(e[e.length-1]));case vt:if(t.target&&this.isFocusableElement(t.target)){const e=t.target;e instanceof qo&&e.childItemLength()>0&&e.expanded?e.expanded=!1:e instanceof qo&&e.parentElement instanceof qo&&qo.focusItem(e.parentElement)}return!1;case gt:if(t.target&&this.isFocusableElement(t.target)){const e=t.target;e instanceof qo&&e.childItemLength()>0&&!e.expanded?e.expanded=!0:e instanceof qo&&e.childItemLength()>0&&this.focusNextNode(1,t.target)}return;case mt:return void(t.target&&this.isFocusableElement(t.target)&&this.focusNextNode(1,t.target));case bt:return void(t.target&&this.isFocusableElement(t.target)&&this.focusNextNode(-1,t.target));case ft:return void this.handleClick(t)}return!0},this.handleSelectedChange=t=>{if(t.defaultPrevented)return;if(!(t.target instanceof Element&&Bo(t.target)))return!0;const e=t.target;e.selected?(this.currentSelected&&this.currentSelected!==e&&(this.currentSelected.selected=!1),this.currentSelected=e):e.selected||this.currentSelected!==e||(this.currentSelected=null)},this.setItems=()=>{const t=this.treeView.querySelector("[aria-selected='true']");this.currentSelected=t,null!==this.currentFocused&&this.contains(this.currentFocused)||(this.currentFocused=this.getValidFocusableItem()),this.nested=this.checkForNestedItems(),this.getVisibleNodes().forEach((t=>{Bo(t)&&(t.nested=this.nested)}))},this.isFocusableElement=t=>Bo(t),this.isSelectedElement=t=>t.selected}slottedTreeItemsChanged(){this.$fastController.isConnected&&this.setItems()}connectedCallback(){super.connectedCallback(),this.setAttribute("tabindex","0"),s.DOM.queueUpdate((()=>{this.setItems()}))}handleClick(t){if(t.defaultPrevented)return;if(!(t.target instanceof Element&&Bo(t.target)))return!0;const e=t.target;e.disabled||(e.selected=!e.selected)}focusNextNode(t,e){const i=this.getVisibleNodes();if(!i)return;const s=i[i.indexOf(e)+t];Be(s)&&qo.focusItem(s)}getValidFocusableItem(){const t=this.getVisibleNodes();let e=t.findIndex(this.isSelectedElement);return-1===e&&(e=t.findIndex(this.isFocusableElement)),-1!==e?t[e]:null}checkForNestedItems(){return this.slottedTreeItems.some((t=>Bo(t)&&t.querySelector("[role='treeitem']")))}getVisibleNodes(){return function(t,e){if(t&&Be(t))return Array.from(t.querySelectorAll(e)).filter((t=>null!==t.offsetParent))}(this,"[role='treeitem']")||[]}}(0,d.gn)([(0,s.attr)({attribute:"render-collapsed-nodes"})],jo.prototype,"renderCollapsedNodes",void 0),(0,d.gn)([s.observable],jo.prototype,"currentSelected",void 0),(0,d.gn)([s.observable],jo.prototype,"slottedTreeItems",void 0);class _o{constructor(t){this.listenerCache=new WeakMap,this.query=t}bind(t){const{query:e}=this,i=this.constructListener(t);i.bind(e)(),e.addListener(i),this.listenerCache.set(t,i)}unbind(t){const e=this.listenerCache.get(t);e&&(this.query.removeListener(e),this.listenerCache.delete(t))}}class Ko extends _o{constructor(t,e){super(t),this.styles=e}static with(t){return e=>new Ko(t,e)}constructListener(t){let e=!1;const i=this.styles;return function(){const{matches:s}=this;s&&!e?(t.$fastController.addStyles(i),e=s):!s&&e&&(t.$fastController.removeStyles(i),e=s)}}unbind(t){super.unbind(t),t.$fastController.removeStyles(this.styles)}}const Wo=Ko.with(window.matchMedia("(forced-colors)")),Go=Ko.with(window.matchMedia("(prefers-color-scheme: dark)")),Yo=Ko.with(window.matchMedia("(prefers-color-scheme: light)"));class Xo{constructor(t,e,i){this.propertyName=t,this.value=e,this.styles=i}bind(t){s.Observable.getNotifier(t).subscribe(this,this.propertyName),this.handleChange(t,this.propertyName)}unbind(t){s.Observable.getNotifier(t).unsubscribe(this,this.propertyName),t.$fastController.removeStyles(this.styles)}handleChange(t,e){t[e]===this.value?t.$fastController.addStyles(this.styles):t.$fastController.removeStyles(this.styles)}}const Qo="not-allowed",Zo=":host([hidden]){display:none}";function Jo(t){return`${Zo}:host{display:${t}}`}const tn=function(){if("boolean"==typeof ze)return ze;if("undefined"==typeof window||!window.document||!window.document.createElement)return ze=!1,ze;const t=document.createElement("style"),e=function(){const t=document.querySelector('meta[property="csp-nonce"]');return t?t.getAttribute("content"):null}();null!==e&&t.setAttribute("nonce",e),document.head.appendChild(t);try{t.sheet.insertRule("foo:focus-visible {color:inherit}",0),ze=!0}catch(t){ze=!1}finally{document.head.removeChild(t)}return ze}()?"focus-visible":"focus"}}]);