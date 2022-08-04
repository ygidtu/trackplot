<template>
  <div>
    <el-row :gutter="20"> <!--  v-if="status.region === 'success'" -->
      <el-col :span="8" :offset="1">
        <el-button @click="dialog.reference = true">Choose reference file</el-button>
      </el-col>
      <el-col :span="12" :offset="1">
        <Param func="set_reference" :path.sync="options.reference"/>
      </el-col>
    </el-row>
    <div id="dialog">
      <el-dialog title="Reference" :visible.sync="dialog.reference" :modal="true">
        <el-row>
          <el-col :span="16" :offset="2">
            <el-input type="textarea"
                      v-model="options.reference"
                      clearable @input="fill_path(options.reference)"
                      :rows="5"
            />
          </el-col>
          <el-col :span="4">
            <el-button type="primary" @click="valid(options.reference)">Choose</el-button>
          </el-col>
        </el-row>

        <el-row>
          <ul class="infinite-list" style="overflow:auto">
            <li v-for="i in options.references" :key="i.path" style="text-align: left;">
              <el-link @click="fill_path(i.path)" :icon="i.isdir ? 'el-icon-folder' : 'el-icon-files'">
                {{ i.path }}
              </el-link>
            </li>
          </ul>
        </el-row>
      </el-dialog>
    </div>
  </div>
</template>

<script>

import Param from '@/components/Param'

export default {
  name: 'Reference',
  components: {
    Param
  },
  data() {
    return {
      dialog: {
        reference: false
      },
      options: {
        references: [],
        reference: ""
      }
    }
  },
  methods: {
    fill_path: function (path) {
      const self = this;

      this.options.reference = path;

      this.axios.get("http://127.0.0.1:5000/api/file", {
        params: {"target": path}
      }).then(response => {
        self.options.references = response.data;
      }).catch(error => {
        self.$notify({
          showClose: true,
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: error.response.data.detail
        });
      })
    },
    valid: function (path) {
      const self = this;
      this.axios.get("http://127.0.0.1:5000/api/file", {
        params: {"target": path, valid: true},
      }).then(response => {
        if (response.data) {
          self.dialog.reference = false;
        } else {
          self.$notify({
            showClose: true,
            type: 'error',
            title: "Error",
            message: "Please select a file, instead of directory"
          });
        }
      }).catch(error => {
        self.$notify({
          showClose: true,
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: error.response.data.detail
        });
      })
    },
  },
  mounted() {
    this.fill_path(this.options.reference)
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
