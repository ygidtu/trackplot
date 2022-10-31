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
        <Param :func.sync="'add_' + this.image.select.toLowerCase()" :path.sync="options.file" :plot_type.sync="image.select"/>
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

<script>
import Param from './Param.vue'
import urls from '../url.js'

import {ref} from 'vue'

export default {
  name: "Add",
  components: {Param},
  data() {
    return {
      image: {
        type: ["Density", "Line", "Heatmap", "IGV", "HiC", "Motif"],
        select: "Density"
      },
      dialog: {
        file: ref(false)
      },
      options: {
        files: [],
        file: ""
      },
    }
  },
  methods: {
    fill_path: function (path) {
      const self = this;

      this.options.file = path;

      this.axios.get(urls.file, {
        params: {"target": path}
      }).then(response => {
        self.options.files = response.data;
      }).catch(error => {
        ElNotification({
          showClose: true,
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: h('i', { style: 'color: teal' },error.response.data.detail)
        })
      })
    },
    valid: function (path) {
      const self = this;
      this.axios.get(urls.file, {
        params: {"target": path, valid: true},
      }).then(response => {
        if (response.data) {
          self.dialog.file = ref(false);
        } else {
         ElNotification({
          showClose: true,
          type: 'error',
          title: `Error`,
          message: h('i', { style: 'color: teal' }, "Please select a file, instead of directory")
        })
        }
      }).catch(error => {
         ElNotification({
          showClose: true,
          type: 'error',
          title: `Error Status: ${error.response.status}`,
          message: h('i', { style: 'color: teal' }, error.response.data.detail)
        })
      })
    }

  },
  mounted() {
    this.fill_path(this.options.file)
  },
}
</script>

<style scoped>

</style>