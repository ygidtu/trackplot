<template>
  <div>
    <el-row :gutter="20">
      <el-col :span="20" :offset="2">
        <el-row>
          <el-col :span="24">
            <el-select v-model="image.select" placeholder="Please choose">
            <el-option v-for="item in image.type" :key="item" :label="item" :value="item"/>
          </el-select>
          </el-col>
        </el-row>
        <el-divider/>
        <param-comp :func.sync="'add_' + image.select.toLowerCase()"
                    :path.sync="options.file"
                    :plot_type.sync="image.select"
                    @select-data="valid" />
      </el-col>
    </el-row>
  </div>
</template>

<script lang="ts" setup>
import axios from 'axios'
import ParamComp from './Param.vue'

interface Path {
  path: string,
  isdir: Boolean
}

interface Option {
  files: Path[],
  file: string
}
</script>

<script lang="ts">
import {errorPrint, Notification} from "../error";

import urls from '../url';
import {h} from 'vue';

export default {
  name: "addComp",
  data() {
    let options: Option = {files: [], file: ""}

    return {
      image: {
        type: ["Density", "Line", "Heatmap", "IGV", "HiC", "Links", "Sites", "Stroke", "Focus"],
        select: "Density"
      },
      options: options
    }
  },
  emits: ["select-data"],
  methods: {
    fill_path (path: string) {
      this.options.file = path;

      axios.get(urls.file, {
        params: {"target": path}
      }).then((response: AxiosResponse) => {
        this.options.files = response.data;
      }).catch((error: AxiosError) => {
        errorPrint(error)
      })
    },
    valid (data: any) {
      axios.get(urls.file, {
        params: {"target": data.path, valid: true},
      }).then((response: AxiosResponse) => {
        if (response.data) {
          data.type = `add_${this.image.select}`
          this.$emit("select-data", data)
        } else {
          let msg: Notification = {
              type: 'error',
              title: `Error`,
              message: h('i', { style: 'color: teal' }, "Please select a file, instead of directory")
          }
          errorPrint(msg)
        }
      }).catch((error: AxiosError) => {
         errorPrint(error)
      })
    },
  },
  mounted() {
    this.fill_path(this.options.file)
  }
}
</script>

<style scoped>

</style>