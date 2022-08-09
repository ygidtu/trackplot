<template>
  <div>
    <h1>{{ msg }}</h1>

    <el-divider/>

    <el-row>
      <el-col :span="16" :offset="2">
        <el-steps align-center
                  finish-status="success"
                  process-status="process"
        >
          <el-step status="wait" title="Set target region"/>
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
      <el-col :span="20" :offset="2">
        <el-form :model="ruleForm" ref="ruleForm" style="width: 100%" label-width="80px" :rules="rules">
          <el-collapse v-model="active">
            <el-collapse-item title="Region" name="0">
              <el-row :gutter="20">
                <el-col :span="20" :offset="2">
                  <el-form-item prop="region">
                    <el-col :span="20">
                      <el-input
                          v-model="ruleForm.region"
                          placeholder="Please input target region, eg: chr1:100-200:+" clearable
                      />
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
                <Reference />
              </div>
            </el-collapse-item>
            <el-collapse-item title="Add" name="2">
              <Add/>
            </el-collapse-item>
            <el-collapse-item title="Draw" name="3">
              <Param func="plot" path="plot"/>
            </el-collapse-item>
          </el-collapse>
        </el-form>
      </el-col>
    </el-row>
  </div>
</template>

<script>

import Add from '../components/Add.vue'
import Param from '../components/Param.vue'
import Reference from '../components/Reference.vue'
import urls from '../url.js'

export default {
  name: "Plot",
  data() {
    let validRegion = (rule, value, callback) => {
      let pattern = /\w+:\d+-\d+:[+-]/i;
      if (!value) {
        return callback(new Error("The region should not be empty!"));
      } else if (!pattern.test(value)) {
        return callback(new Error("The input region format is wrong!"));
      }
      callback();
    };
    return {
      msg: "Make your own plot",
      active: "0",
      dialog: {
        reference: false
      },
      image: ["Density", "Line", "Heatmap", "IGV"],
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
      plot: null
    };
  },
  methods: {
    remoteFilePath(query) {
      if (query !== "") {
        this.loading = true;
        setTimeout(() => {
          this.loading = false;
          this.options = this.list.filter(item => {
            return item.label.toLowerCase()
                .indexOf(query.toLowerCase()) > -1;
          });
        }, 200);
      } else {
        this.options = [];
      }
    },
    submitRegion() {
      let regions = this.ruleForm.region.split(":")
      let chrom = regions[0]
      let strand = regions[regions.length - 1]
      let sites = regions[1].split("-")

      const self = this;
      this.axios.post(`${urls.plot}?pid=${this.$cookies.get("plot")}&func=set_region`, {
        path: "",
        param: [
          {key: "chromosome", default: chrom, annotation: "str"},
          {key: "start", default: parseInt(sites[0], 10), annotation: "int"},
          {key: "end", default: parseInt(sites[1], 10), annotation: "int"},
          {key: "strand", default: strand, annotation: "str"},
        ]
      }).then(response => {
        ElNotification({
          title: 'Success',
          message: `set_region execute success`,
          type: 'success'
        })
      }).catch(error => {
        ElNotification({
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: h('i', { style: 'color: teal' }, error.response.data.detail)
        })
      })
    },
    reset() {
      this.axios.get(`${urls.del}?pid=${this.$cookies.get("plot")}`)
      location.reload()
    }
  },
  components: {Add, Param, Reference},
  mounted() {
    if (this.$cookies.isKey("plot")) {
      this.axios.get(`${urls.del}?pid=${this.$cookies.get("plot")}`)
    }
    this.$cookies.set("plot", (Math.random() + 1).toString(36).substring(7)) //
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
