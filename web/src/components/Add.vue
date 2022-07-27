<template>
  <div>
    <el-row :gutter="20"> <!--  v-if="status.region === 'success'" -->
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
        <Param :func.sync="funcName" :path.sync="options.file" :plot_type.sync="image.select"/>
      </el-col>
    </el-row>
    <div id="dialog">
      <el-dialog title="Reference" :visible.sync="dialog.file" :modal="true">
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

      </el-dialog>
    </div>
  </div>
</template>

<script>
import Param from '@/components/Param'

export default {
  name: "Add",
  components: {Param},
  data() {
    return {
      image: {
        type: ["Density", "Line", "Heatmap", "IGV"],
        select: "Density"
      },
      dialog: {
        file: false
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

      this.axios.get("/api/file", {
        params: {"target": path}
      }).then(response => {
        self.options.files = response.data;
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
      this.axios.get("/api/file", {
        params: {"target": path, valid: true},
      }).then(response => {
        if (response.data) {
          self.dialog.file = false;
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
    }

  },
  mounted() {
    this.fill_path(this.options.file)
  },
  computed: {
    funcName() {
      return "add_" + this.image.select.toLowerCase()
    }
  }
}
</script>

<style scoped>

</style>