import Vue from 'vue'

import axios from 'axios';
import VueAxios from 'vue-axios';
import VueCookies from 'vue-cookies'

import ElementUI from 'element-ui';
import 'element-ui/lib/theme-chalk/index.css';
import App from './App.vue'
import index from './router';


Vue.use(ElementUI);
Vue.use(VueAxios, axios);
Vue.config.productionTip = false
Vue.use(VueCookies, { expire: '7d'})


new Vue({
  router: index,
  render: h => h(App),
}).$mount('#app')
