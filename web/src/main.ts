import { createApp} from 'vue'
import App from './App.vue'
import { createRouter, createWebHashHistory, RouteRecordRaw } from 'vue-router';
import { VueCookieNext } from 'vue-cookie-next'
import VueAxios from 'vue-axios'
import axios from 'axios'
import 'element-plus/dist/index.css'


const routes: RouteRecordRaw[] = [
    { path: "/", component: () => import("./pages/Home.vue") },
    { path: "/home", redirect: "/" },
    { path: "/plot", component: () => import("./pages/Plot.vue") },
    // { path: "/log", component: () => import("./components/Img.vue") },
];

const router = createRouter({
    history: createWebHashHistory(),
    routes,
})

// set default config
VueCookieNext.config({ expire: '7d' })

// set global cookie
VueCookieNext.setCookie('theme', 'default')

const app = createApp(App)
app.use(router)
app.use(VueCookieNext)
app.use(VueAxios, axios)
app.provide('axios', app.config.globalProperties.axios)  // provide 'axios'
app.mount('#app');
