<script lang="ts" setup>
import {Delete, View} from "@element-plus/icons-vue"
import AddComp from '../components/Add.vue'
import Annotation from '../components/Annotation.vue'
import ParamComp from '../components/Param.vue'
import LogComp from "../components/Log.vue"
</script>

<template>
  <div>
    <h1>{{ msg }}</h1>

    <el-divider/>

    <el-row>
      <el-col :span="16" :offset="2">
        <el-steps align-center
                  finish-status="success"
                  process-status="process"
                  :active="active">
          <el-step title="Set target region"/>
          <el-step title="Set annotation"/>
          <el-step title="Set plot details"/>
          <el-step title="Preview/Save"/>
        </el-steps>
      </el-col>
      <el-col :span="4" :offset="2">
        <el-button type="danger" :icon="Delete" @click="reset">Reset</el-button>
      </el-col>
    </el-row>
    <el-divider/>
    <el-row :gutter="20" style="overflow: hidden">
      <el-col :span="12">
        <el-form :model="ruleForm" ref="ruleForm" style="width: 100%" label-width="80px" :rules="rules">
          <el-tabs v-model="active">
            <el-tab-pane label="Region" :name="0">
              <el-row :gutter="20">
                <el-col :span="20" :offset="2">
                  <el-row :gutter="10">
                    <el-col :span="20">
                      <el-input
                        v-model="ruleForm.region" clearable
                        placeholder="Please input target region, eg: chr1:100-200:+" />
                    </el-col>
                    <el-col :span="4">
                      <el-button type="primary" @click="submitRegion">Confirm</el-button>
                    </el-col>
                  </el-row>
                </el-col>
              </el-row>
            </el-tab-pane>
            <el-tab-pane label="Annotation" :name="1">
              <el-scrollbar :max-height="windowHeight" always>
                <el-col :span="24">
                  <annotation @select-data="makeProgress"/>
                </el-col>
              </el-scrollbar>
            </el-tab-pane>
            <el-tab-pane label="Add" :name="2" always>
              <el-scrollbar :max-height="windowHeight">
                <el-col :span="24">
                  <add-comp  @select-data="makeProgress" />
                </el-col>
              </el-scrollbar>
            </el-tab-pane>
            <el-tab-pane label="Draw" :name="3" always>
              <el-scrollbar :max-height="windowHeight">
                <el-col :span="24">
                  <param-comp func="plot" path="plot" @select-data="makeProgress" />
                </el-col>
              </el-scrollbar>
            </el-tab-pane>
          </el-tabs>
        </el-form>
      </el-col>
      <el-col :span="12">
        <el-tabs v-model="setting">
          <el-tab-pane label="Plot settings" name="settings">
            <el-scrollbar :height="windowHeight">
              <el-space direction="vertical" fill>
                <el-row v-if="progress.region !== null">
                  <el-descriptions class="margin-top" title="Plotting region" :column="4" border>
                    <el-descriptions-item v-for="p in progress.region.param" :label="p.key" :key="p.key">
                      {{ p.default }}
                    </el-descriptions-item>
                  </el-descriptions>
                </el-row>

                <el-row v-if="progress.annotation !== null" >
                  <el-col :span="24">
                    <el-descriptions class="margin-top" title="Reference" :column="4" border>
                    <el-descriptions-item label="Path">
                      <el-popover
                        placement="top-start"
                        title="Detailed"
                        :width="'70%'"
                        trigger="hover"
                      >
                        <template #reference>
                          <el-text type="info">{{ progress.annotation.path }}<el-icon><View /></el-icon></el-text>
                        </template>
                        <el-descriptions title="Parameters" :column="3" border>
                          <el-descriptions-item
                              v-for="p in showParam(progress.annotation.param)"
                              :key="p.key" :label="p.key">{{ p.default }}
                          </el-descriptions-item>
                        </el-descriptions>
                      </el-popover>
                    </el-descriptions-item>
                  </el-descriptions>
                  </el-col>
                </el-row>

                <el-row v-if="progress.files.length > 0" >
                  <el-col :span="24">
                    <el-descriptions class="margin-top" title="Files" :column="2" border>
                      <div v-for="(file, idx) in progress.files" :key="`${file.path}:${file.type}`">
                          <el-descriptions-item :label="file.type.replace(/add_/, '')">
                            <el-popover
                              placement="top-start"
                              title="Detailed"
                              :width="'70%'"
                              trigger="hover"
                            >
                              <template #reference><el-text type="info">{{ file.path }}<el-icon><View /></el-icon></el-text></template>
                              <el-descriptions title="Parameters" :column="3" border>
                                <el-descriptions-item v-for="p in showParam(file.param)" :key="p.key" :label="p.key">
                                  {{ p.default }}
                                </el-descriptions-item>
                              </el-descriptions>
                            </el-popover>
                        </el-descriptions-item>
                        <el-descriptions-item label="Remove this">
                          <el-button type="danger" :icon="Delete" circle @click="removeThis(idx)" />
                        </el-descriptions-item>
                      </div>
                    </el-descriptions>
                  </el-col>
                </el-row>

                <el-row v-if="progress.draw !== null" >
                  <el-col :span="24">
                    <el-descriptions class="margin-top" title="Plotting parameters" :column="4" border>
                      <el-descriptions-item :label="progress.draw.type">
                      <el-popover
                          placement="top-start"
                          title="Detailed"
                          :width="'70%'"
                          trigger="hover"
                        >
                          <template #reference><el-text type="info">View<el-icon><View /></el-icon></el-text></template>
                          <el-descriptions title="Parameters" :column="3" border>
                            <el-descriptions-item v-for="p in showParam(progress.draw.param)" :key="p.key" :label="p.key">
                              {{ p.default }}
                            </el-descriptions-item>
                          </el-descriptions>
                        </el-popover>
                      </el-descriptions-item>
                    </el-descriptions>
                  </el-col>
                </el-row>
              </el-space>
            </el-scrollbar>
          </el-tab-pane>
          <el-tab-pane label="Process log" name="log" v-if="progress.draw !== null">
            <el-scrollbar :height="windowHeight">
              <el-col :span="24">
                <log-comp v-if="pid !== ''" :pid="pid" />
              </el-col>
            </el-scrollbar>
          </el-tab-pane>
        </el-tabs>
      </el-col>
    </el-row>


    <el-dialog
      v-model="showImage"
      title="Preview"
      width="60%"
      :before-close="handleClose"
    >
      <template #footer>
        <el-row>
          <el-col :span="20" :offset="2">
            <el-image :src="img"></el-image>
          </el-col>
        </el-row>
      </template>
    </el-dialog>
  </div>
</template>

<script lang="ts">

import {defineComponent} from "vue";
import {AxiosRequestConfig, AxiosResponse, AxiosError} from "axios";
import {saveAs} from "file-saver";
import urls from '../url'
import {errorPrint, Notification} from "../error";
import {ElLoading} from "element-plus";


interface Param {
  key: string,
  annotation: string,
  default: any,
  note: string | null
}

interface FilePath {
  path: string,
  type: string,
  param: Param[]
}

interface Progress {
  region: string | null,
  annotation: FilePath | null,
  files: FilePath[] | null,
  draw: FilePath | null
}

const validRegion = (_: any, value: any, callback: any) => {
  let pattern = /\w+:\d+-\d+:[+-]/i;
  if (!value) {
    return callback(new Error("The region should not be empty!"));
  } else if (!pattern.test(value)) {
    return callback(new Error("The input region format is wrong!"));
  }
  callback();
}

export default defineComponent({
  name: "Plot",
  data() {
    let progress: Progress = {
      region: null,
      annotation: null,
      files: [],
      draw: null
    }
    let img: Blob | null | string = null
    return {
      msg: "Make your own plot",
      active: 0, setting: "settings", pid: "",
      ruleForm: { region: "chr1:1270656-1284730:+" },
      rules: {
        region: [
          {validator: validRegion, trigger: "blur"}
        ]
      },
      progress: progress,
      img: img,
    };
  },
  methods: {
    submitRegion() {
      let regions = this.ruleForm.region.split(":")
      let chrom = regions[0]
      let strand = regions[regions.length - 1]
      let sites = regions[1].split("-")

      this.progress.region = {
        path: "", type: "set_region",
        param: [
            {key: "chromosome", default: chrom, annotation: "str"},
            {key: "start", default: parseInt(sites[0], 10), annotation: "int"},
            {key: "end", default: parseInt(sites[1], 10), annotation: "int"},
            {key: "strand", default: strand, annotation: "str"}
        ]
      }
    },
    reset() {
      this.axios.get(`${urls.del}?pid=${this.$cookie.getCookie("plot")}`)
      location.reload()
    },
    showParam(params: Param[]) {
      let p = []

      for (let param of params) {
        if (param.default !== null && param.default.toString().match(/empty/) === null) {
          p.push(param)
        }
      }
      return p
    },
    handleClose() { this.img = null },
    makeProgress(file: FilePath) {
      if (file.type === "annotation") {
        this.progress.annotation = file
      } else if (file.type !== "plot" && file.type !== "save") {
        this.progress.files?.push(file)
      } else if (file.type === "plot" || file.type == "save") {
        this.progress.draw = file
        this.setting = "log"
        this.submit(file.type)
      }
    },
    submit(type: string) {
      const loading = ElLoading.service({
        lock: true,
        text: 'Loading',
        background: 'rgba(0, 0, 0, 0.7)',
      })

      let config: AxiosRequestConfig = {responseType: "application/json"}
      if (type === "plot" || type === "save") {
        config["responseType"] = "blob"
        this.img = null
      } else if (type === "") {
        let msg: Notification = {
          type: 'error',
          title: `Please set up correct file path`,
          message: ""
        }
        errorPrint(msg)
        return
      }

      this.axios.post(
          `${urls.plot}/${this.$cookie.getCookie("plot")}`,
          this.progress, config
      ).then((response: AxiosResponse) => {
        if (type === "plot") {
          const blob = new Blob([response.data], {type: response.headers['content-type']})
          this.img = window.URL.createObjectURL(blob)
        } else if (type === "save") {
          let filename = response.headers["content-disposition"].split("filename=")[1]
          saveAs(response.data, filename)
        } else {
          let msg: Notification = {
            title: 'Success',
            message: `${type} execute success`,
            type: 'success'
          }
          errorPrint(msg)
        }
      }).catch((error: AxiosError) => {
        errorPrint(error)
      }).finally(() => {
        loading.close()
      })
    },
    removeThis(idx: Number) {
      let kept = [];
      for (let i = 0; i < this.progress.files.length; i++) {
        if (i !== idx) {
          kept.push(this.progress.files[i])
        }
      }
      this.progress.files = kept
    }
  },
  computed: {
    showImage() { return this.img !== null },
    windowHeight: function() {
      return window.innerHeight - 300
    }
  },
  mounted() {
    if (this.$cookie.isCookieAvailable("plot")) {
      this.axios.get(`${urls.del}?pid=${this.$cookie.getCookie("plot")}`)
    }

    this.$cookie.setCookie("plot", (Math.random() + 1).toString(36).substring(7))
    // this.$cookie.setCookie("plot", "test")
    this.pid = this.$cookie.getCookie("plot")
  },
  beforeUnmount() {
      this.axios.get(`${urls.del}?pid=${this.$cookie.getCookie("plot")}`)
  },
})
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
