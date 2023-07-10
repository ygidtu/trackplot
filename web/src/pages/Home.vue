<template>
  <div>
    <h1>{{ msg }}</h1>

    <el-divider/>

    <el-row>
      <el-col :span="20" :offset="2">
        <el-link type="primary" href="/#/plot">Create your own plot</el-link>
      </el-col>
    </el-row>

    <el-divider/>

    <el-row>
      <el-col :span="20" :offset="2">
        <el-image :src="example" width="100%"/>
      </el-col>
    </el-row>
  </div>
</template>

<script lang="ts">

export default {
  name: 'Home',
  data() {
    return {
      msg: "Welcome to trackplot App",
      example: `https://trackplot.readthedocs.io/en/latest/imgs/diagram.png`
    }
  }
}
</script>

<!-- Add "scoped" attribute to limit CSS to this component only -->
<style scoped>
h3 {
  margin: 40px 0 0;
}

ul {
  list-style-type: none;
  padding: 0;
}

li {
  display: inline-block;
  margin: 0 10px;
}

a {
  color: #42b983;
}
</style>
