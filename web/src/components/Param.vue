<template>
  <div class="params">
    <el-form ref="form" v-model="param" label-width="120px">
      <div v-for="(p, index) in param" :key="index">
        <el-form-item :label="p.key" v-if="!String(p.default).includes('inspect._empty')">
          <el-input v-if="p.annotation === 'str' || p.annotation === 'Optional[str]'" v-model="p.default"/>
          <el-input-number v-else-if="p.annotation === 'int'" v-model="p.default"/>
          <el-input-number v-else-if="p.annotation === 'float'" :precision='2' v-model="p.default"/>
          <el-switch v-else-if="p.annotation === 'bool'" v-model="p.default" active-text="True" inactive-text="False"/>
          <el-input v-else-if="!p.default.includes('inspect._empty')" v-model="p.default"></el-input>
        </el-form-item>
      </div>

      <el-button type="primary" @click="submit" v-if="func !== 'plot'">Confirm</el-button>
      <el-button-group v-else>
        <el-button type="primary" icon="el-icon-view" @click="submit">Preview</el-button>
        <el-button type="primary" @click="save">Save pdf<i class="el-icon-download"></i></el-button>
      </el-button-group>
    </el-form>

    <div v-if="img !== null">
      <el-divider/>
      <el-row>
        <el-col :span="20" :offset="2">
          <el-image :src="img"></el-image>
        </el-col>
      </el-row>
    </div>
  </div>
</template>

<script>
import saveAs from 'file-saver';
import urls from '../url.js';

export default {
  name: "Param",
  props: {
    func: {required: true, type: String},
    path: {required: true, type: String},
    plot_type: {type: String}
  },
  data() {
    return {
      param: {},
      img: null
    }
  },
  methods: {
    loadParams() {
      const self = this;
      this.axios.get(urls.params, {
        params: {target: this.$props.func}
      }).then(response => {
        self.param = []
        for (let row of response.data) {
          if (row.annotation === 'float') {
            row.default = parseFloat(row.default)
          } else if (row.annotation === 'int') {
            row.default = parseInt(row.default)
          }

          self.param.push(row)
        }
      }).catch(error => {
          ElNotification({
            title: 'Title',
            message: h('i', { style: 'color: teal' }, 'This is a reminder'),
          })
        self.$notify({
          showClose: true,
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: error.response.data
        });
      })
    },
    submit() {
      const self = this;

      let config = {}
      if (this.$props.func === "plot") {
        config["responseType"] = "blob"
      } else if (this.$props.path === "") {
        ElNotification({
           showClose: true,
           type: 'error',
           title: `Please choose file first`,
        })
        return
      }

      this.axios.post(`${urls.plot}?pid=${this.$cookies.get("plot")}&func=${this.$props.func}`, {
        path: this.$props.path,
        param: this.param
      }, config).then(response => {
        if (this.$props.func === "plot") {
          const {data, headers} = response
          const blob = new Blob([data], {type: headers['content-type']})
          self.img = window.URL.createObjectURL(blob)
        } else {
          ElNotification({
            title: 'Success',
            message: `${this.$props.func} execute success`,
            type: 'success'
          })
        }
      }).catch(error => {
        ElNotification({
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: error.response.data.detail
        })
      })
    },
    save() {
      const self = this;
      if (this.$props.func !== "plot") {
        ElNotification({
          type: 'error',
          title: `Please choose file first`,
        })
        return
      }

      this.axios.post(`${urls.plot}?pid=${this.$cookies.get("plot")}&func=save`, {
        path: this.$props.path,
        param: this.param,
      }, {responseType: 'blob'}).then(response => {
        const {data, headers} = response
        let filename = headers["content-disposition"]
        saveAs(data, filename)
      }).catch(error => {
        ElNotification({
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: error.response.data.detail
        })
      })
    },
  },
  mounted() {
    this.loadParams()
  },
  watch: {
    func: {
      handler: function () {
        this.loadParams()
      },
      deep: true
    }
  }
}
</script>

<style scoped>

</style>