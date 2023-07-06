<template>
  <div>
    <el-row :gutter="20">
      <el-col :span="8" :offset="1">
        <el-row>
          <el-select v-model="image.select" placeholder="Please choose">
            <el-option v-for="item in image.type" :key="item" :label="item" :value="item">
            </el-option>
          </el-select>
        </el-row>
        <el-divider/>
        <el-row>
          <el-button @click="dialog.file = true">Choose file</el-button>
        </el-row>
      </el-col>
      <el-col :span="12" :offset="1">
        <param-comp :func.sync="'add_' + image.select.toLowerCase()"
                    :path.sync="options.file"
                    :plot_type.sync="image.select"
                    @select-data="selectData"
        />
      </el-col>
    </el-row>
    <div id="dialog">
      <el-dialog title="Reference" v-model="dialog.file">
        <template #footer>
          <el-row>
            <el-col :span="16" :offset="2">
              <el-input type="textarea"
                        v-model="options.file"
                        clearable @input="fill_path(options.file)"
                        :rows="5"
              />
            </el-col>
            <el-col :span="4">
              <el-button type="primary" @click="valid(options.file)">Choose</el-button>
            </el-col>
          </el-row>

          <el-row>
            <ul class="infinite-list" style="overflow:auto">
              <li v-for="i in options.files" :key="i.path" style="text-align: left;">
                <el-link @click="fill_path(i.path)" :icon="i.isdir ? 'el-icon-folder' : 'el-icon-files'">
                  {{ i.path }}
                </el-link>
              </li>
            </ul>
          </el-row>
        </template>
      </el-dialog>
    </div>
  </div>
</template>

<script lang="ts" setup>
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
import {errorPrint, Notification} from "../error.ts";
import {AxiosError, AxiosResponse} from "axios";
import urls from '../url';
import {h} from 'vue';

export default {
  name: "addComp",
  data() {
    let options: Option = {files: [], file: ""}

    return {
      image: {
        type: ["Density", "Line", "Heatmap", "IGV", "HiC", "Motif"],
        select: "Density"
      },
      dialog: {file: false},
      options: options
    }
  },
  emits: ["select-data"],
  methods: {
    fill_path (path: string) {
      this.options.file = path;

      this.axios.get(urls.file, {
        params: {"target": path}
      }).then((response: AxiosResponse) => {
        this.options.files = response.data;
      }).catch((error: AxiosError) => {
        errorPrint(error)
      })
    },
    valid (path: string) {
      this.axios.get(urls.file, {
        params: {"target": path, valid: true},
      }).then((response: AxiosResponse) => {
        if (response.data) {
          this.dialog.file = false;
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
    selectData(data: any) {
      data.type = this.image.type
      this.$emit("select-data", data)
    }
  },
  mounted() {
    this.fill_path(this.options.file)
  }
}
</script>

<style scoped>

</style>