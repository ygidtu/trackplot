import { createApp, ref } from 'vue'
import App from './App.vue'
import axios from 'axios'
import VueAxios from 'vue-axios'
import VueCookies from 'vue-cookies'
import { createRouter, createWebHashHistory } from 'vue-router'
import {ElNotification} from "element-plus";
import * as ElementPlusIconsVue from '@element-plus/icons-vue'

import Home from "./pages/Home.vue";
import Plot from "./pages/Plot.vue";


const routes = [
    { path: "/", component: Home },
    { path: "/home", redirect: "/" },
    { path: "/plot", component: Plot },
];

const router = createRouter({
    history: createWebHashHistory(),
    routes,
})

const app = createApp(App);
app.use(ElNotification)
app.use(ref)
app.use(VueAxios, axios)
app.use(VueCookies)
app.use(router)

for (const [key, component] of Object.entries(ElementPlusIconsVue)) {
  app.component(key, component)
}

app.mount('#app')
