<template>
  <el-row>
    <el-col :span="20" :offset="2">
        <router-view style="width:100%" />
    </el-col>
  </el-row>

</template>

<script>

export default {
  name: 'App'
}
</script>

<style>
#app {
  font-family: Avenir, Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
}
</style>
