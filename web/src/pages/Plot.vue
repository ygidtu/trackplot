<script lang="ts" setup>
import AddComp from '../components/Add.vue'
import ParamComp from '../components/Param.vue'
import Reference from '../components/Reference.vue'
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
                  :active="activeSteps"
        >
          <el-step title="Set target region"/>
          <el-step title="Set reference"/>
          <el-step title="Set plot details"/>
          <el-step title="Preview/Save"/>
        </el-steps>
      </el-col>
      <el-col :span="4" :offset="2">
        <el-button type="danger" icon="el-icon-delete" @click="reset">Reset</el-button>
      </el-col>
    </el-row>
    <el-divider/>
    <el-row>
      <el-col :span="12">
        <el-form :model="ruleForm" ref="ruleForm" style="width: 100%" label-width="80px" :rules="rules">
          <el-collapse v-model="active">
            <el-collapse-item title="Region" name="0">
              <el-row :gutter="20">
                <el-col :span="20" :offset="2">
                  <el-form-item prop="region">
                    <el-col :span="20">
                      <el-input
                        v-model="ruleForm.region" clearable
                        placeholder="Please input target region, eg: chr1:100-200:+" />
                    </el-col>
                    <el-col :span="4">
                      <el-button type="primary" @click="submitRegion">Confirm</el-button>
                    </el-col>
                  </el-form-item>
                </el-col>
              </el-row>
            </el-collapse-item>
            <el-collapse-item title="Reference" name="1">
              <div>
                <reference @select-data="makeProgress"/>
              </div>
            </el-collapse-item>
            <el-collapse-item title="Add" name="2">
              <add-comp />
            </el-collapse-item>
            <el-collapse-item title="Draw" name="3">
              <param-comp func="plot" path="plot" @select-data="makeProgress"/>
            </el-collapse-item>
          </el-collapse>
        </el-form>
      </el-col>
      <el-col :span="12">
        <el-divider>Process log</el-divider>
        <log-comp v-if="pid !== ''" :pid="pid" />
      </el-col>
    </el-row>
  </div>
</template>

<script lang="ts">

import {defineComponent, h} from "vue";

import urls from '../url.ts'
import {errorPrint, Notification} from "../error.ts";

interface FilePath {
  path: string,
  type: string,
  param: any
}

interface Progress {
  region: string | null,
  reference: FilePath | null,
  files: FilePath[] | null,
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
      reference: null,
      files: []
    }
    return {
      msg: "Make your own plot",
      pid: "",
      active: "0",
      dialog: {
        reference: false
      },
      image: ["Density", "Line", "Heatmap", "IGV", "HiC", "Motif"],
      options: {
        references: []
      },
      ruleForm: {
        region: "chr1:1270656-1284730:+",
      },
      rules: {
        region: [
          {validator: validRegion, trigger: "blur"}
        ]
      },
      plot: null,
      loading: false,
      progress: progress
    };
  },
  methods: {
    submitRegion() {
      let regions = this.ruleForm.region.split(":")
      let chrom = regions[0]
      let strand = regions[regions.length - 1]
      let sites = regions[1].split("-")

      this.axios.post(`${urls.plot}?pid=${this.$cookie.getCookie("plot")}&func=set_region`, {
        path: "",
        param: [
          {key: "chromosome", default: chrom, annotation: "str"},
          {key: "start", default: parseInt(sites[0], 10), annotation: "int"},
          {key: "end", default: parseInt(sites[1], 10), annotation: "int"},
          {key: "strand", default: strand, annotation: "str"},
        ]
      }).then(() => {
        let msg: Notification = {
          title:'Success',
          message:`set_region execute success`,
          type:'success'
        }
        errorPrint(msg)

        this.progress.region = this.ruleForm.region
      }).catch((error: any) => {
        errorPrint(error)
      })
    },
    reset() {
      this.axios.get(`${urls.del}?pid=${this.$cookie.getCookie("plot")}`)
      location.reload()
    },
    makeProgress(file: FilePath) {
      if (file.type === "reference") {
        this.progress.reference = file
      } else {
        this.progress.files?.push(file)
      }
    },
    submit() {
      const self = this;
      let config: AxiosRequestConfig = {responseType: "json"}
      if (this.$props.func === "plot") {
        config["responseType"] = "blob"
      } else if (this.$props.path === "") {
        let msg: Notification = {
          type: 'error',
          title: `Please set up correct file path`,
          message: ""
        }
        errorPrint(msg)
        return
      }

      this.axios.post(
          `${urls.plot}?pid=${this.$cookie.getCookie("plot")}&func=${this.$props.func}`,
          {
            path: this.$props.path,
            param: this.param
          },
          config
      ).then((response: AxiosResponse) => {
        if (this.$props.func === "plot") {
          const {data, headers} = response
          const blob = new Blob([data], {type: headers['content-type']})
          self.img = window.URL.createObjectURL(blob)
        } else {
          let msg: Notification = {
            title: 'Success',
            message: `${this.$props.func} execute success`,
            type: 'success'
          }
          errorPrint(msg)
        }
      }).catch((error: AxiosError) => {
        errorPrint(error)
      })
    },
  },
  computed: {
    activeSteps() {
      if (this.progress.region === null || this.progress.region === "") {
        return 0
      } else if (this.progress.reference === null) {
        return 1
      } else if (this.progress.files === null || this.progress.files.length < 1) {
        return 2
      } else {
        return 3
      }
    }
  },
  mounted() {
    if (this.$cookie.isCookieAvailable("plot")) {
      this.axios.get(`${urls.del}?pid=${this.$cookie.getCookie("plot")}`)
    }

    // this.$cookie.setCookie("plot", (Math.random() + 1).toString(36).substring(7))
    this.$cookie.setCookie("plot", "test")
    this.pid = this.$cookie.getCookie("plot")
  }
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
