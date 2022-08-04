import Vue from 'vue'
import Router from 'vue-router'
import Home from "@/pages/Home";
import Plot from "@/pages/Plot";

Vue.use(Router);


const routes = [
    { path: "/", component: Home },
    { path: "/home", redirect: "/" },
    { path: "/plot", component: Plot },
];


const originalPush = Router.prototype.push
   Router.prototype.push = function push(location) {
   return originalPush.call(this, location).catch(err => err)
}


const index = new Router({
    routes: routes,
});


window.router = index;

export default index
