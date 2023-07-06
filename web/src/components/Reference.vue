<template>
  <div>
    <el-row :gutter="20">
      <el-col :span="24" :offset="1">
        <param-comp func="set_reference" @select-data="valid" :postfix="/.(gtf|gff\d?)(.gz)?$/" />
      </el-col>
    </el-row>
  </div>
</template>

<script lang="ts" setup>
import ParamComp from './Param.vue'
</script>

<script lang="ts">
import {h} from 'vue'
import urls from '../url';
import {errorPrint, Notification} from "../error";
import {AxiosError} from "axios";

export default {
  name: "reference",
  data() {
    return {
      dialog: {reference: false, path: ""},
    }
  },
  emits: ["select-data"],
  methods: {
    valid (data: any) {
      this.axios.get(urls.file, {
        params: {"target": data.path, valid: true},
      }).then((response: any) => {
        if (response.data) {
          this.dialog.reference = false;
          data.type = "reference"
          this.submit(data)
        } else {
          let msg: Notification = {
            type: 'error',
            title: "Error",
            message: h('i', { style: 'color: teal' }, "Please select a file, instead of directory")
          }
          errorPrint(msg)
        }
      }).catch((error: any) => {
        errorPrint(error)
      })
    },
    submit(data: any) {
      this.axios.post(
          `${urls.plot}?pid=${this.$cookie.getCookie("plot")}&func=set_reference`,
          data,
      ).then((_) => {
        let msg: Notification = {
          title: 'Success',
          message: `set_reference execute success`,
          type: 'success'
        }
        errorPrint(msg)
        this.$emit("select-data", data)
      }).catch((error: AxiosError) => {
        errorPrint(error)
      })
    },
  }
}
</script>

<!-- Add "scoped" attribute to limit CSS to this component only -->
<style scoped>
h3 {
  margin: 40px 0 0;
}

a {
  color: #42b983;
}
</style>
