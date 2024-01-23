<script setup lang="ts">
import {Message, Document, Folder, View, Download} from '@element-plus/icons-vue'
</script>

<template>
  <div class="params">
    <el-form ref="form" label-width="180px">
      <div v-for="(p, index) in param" :key="index">
        <el-row :gutter="20" v-if="String(p.default).match(/empty/) === null">
          <el-col :span="20">
            <el-form-item :label="p.key" >
              <el-input v-if="p.annotation === 'str' || p.annotation === 'Optional[str]'" v-model="p.default"/>
              <el-input-number v-else-if="p.annotation === 'int'" v-model="p.default"/>
              <el-input-number v-else-if="p.annotation === 'float'" :precision='2' v-model="p.default" :step=".1"/>
              <el-color-picker v-else-if="p.annotation === 'color'" v-model="p.default" />

              <el-radio-group v-model="p.default" v-else-if="p.annotation === 'choice'">
                <el-radio v-for="i in p.choice" type="primary" :key="i" :label="i">{{i}}</el-radio>
              </el-radio-group>

              <el-select v-model="p.default" v-else-if="p.annotation === 'select'">
                <el-option v-for="i in p.choice" type="primary" :key="i" :value="i" :label="i">{{i}}</el-option>
              </el-select>

              <el-switch v-else-if="p.annotation === 'bool'" v-model="p.default" active-text="True" inactive-text="False"/>
              <el-input v-else-if="p.default.match(/empty/) === null" v-model="p.default"></el-input>
            </el-form-item>
          </el-col>
          <el-col :span="4">
            <el-popover v-if="p.note !== ''"
                  placement="top-start"
                  title="Description"
                  :width="400"
                  trigger="hover"
                >
              <template #reference>
                <el-button :icon="Message" circle />
              </template>
              <el-text type="info" v-if="!p.note.toString().match(/\n/)">{{ p.note }}</el-text>
              <div v-else>
                <ul>
                  <el-text type="info" v-for="(content, index) in p.note.split(/\n/)">
                    <li v-if="index > 0"> <el-text type="danger">{{content.split(":")[0]}}</el-text>:</li> {{ content.split(":")[1] }}
                  </el-text>
                </ul>
              </div>
            </el-popover>
          </el-col>
        </el-row>
      </div>

      <div v-if="func.match(/(sites|links|stroke|focus|plot)/) === null">
        <el-form-item label="File path">
          <el-input v-model="file" @input="fill_path(file)" placeholder="please select or input the file path" />
        </el-form-item>

        <el-divider />

        <ul class="infinite-list" style="overflow: auto">
          <li v-for="i in files.Dirs" :key="i.path" class="infinite-list-item" @click="fill_path(i.path)">
            <el-link type="primary" :icon="Folder">{{ i.path }}</el-link>
          </li>
          <li v-for="i in files.Files" :key="i.path" class="infinite-list-item" @click="fill_path(i.path)">
            <el-link type="primary" :icon="Document">{{ i.path }}</el-link>
          </li>
        </ul>
      </div>

      <el-button type="primary" @click="submit" v-if="func !== 'plot'">Confirm</el-button>
      <el-button-group v-else>
        <el-button type="primary" :icon="View" @click="submit">Preview</el-button>
        <el-button type="primary" :icon="Download" @click="save">Save pdf</el-button>
      </el-button-group>
    </el-form>
  </div>
</template>


<script lang="ts">
import {defineComponent} from 'vue'
import axios from 'axios'
import {AxiosResponse, AxiosError} from 'axios'
import urls from '../url'

import { errorPrint } from "../error.ts";


interface Param {
  key: String,
  annotation: String,
  default: String,
  note: String,
  choice: string
}

interface Path {
  path: string,
  isdir: Boolean
}

interface Files {
  Dirs: Path[],
  Files: Path[]
}

interface Data {
  param: Param[],
  file: string,
  files: Files
}

const decodeChoice = (choices: String) => {
  let res = Array<String>();
  let choiceMatch = choices.match(/(choice|select)\[(.*)\]/);
  if (choiceMatch === null) {
    return res
  }
  
  let choice = choiceMatch[2].toString().replaceAll("'", "").replaceAll('"', '').split(',');

  for (let c of choice) {
    res.push(c.trim())
  }
  return res
}

export default defineComponent({
  props: {
    func: {required: true, type: String},
    plot_type: {type: String},
    postfix: {type: RegExp, default: null}
  },
  data() {
    let d: Data = {
      param: [], file: "",
      files: {Dirs: [], Files: []}
    }
    return d
  },
  emits: ["select-data"],
  methods: {
    save() {
      this.$emit("select-data", {path: this.file, param: this.param, type: "save"})
    },
    preview() {
      this.file = "preview"
      this.submit()
    },
    submit() {
      this.$emit("select-data", {path: this.file, param: this.param, type: this.$props.func})
    },
    loadParams() {
      // 👇️ const data: GetUsersResponse
      axios.get(urls.params, {params: {target: this.$props.func}}
      ).then((response: AxiosResponse) => {
        this.param = []
        for (let row of response.data) {
          if (row.annotation === 'float') {
            row.default = parseFloat(row.default)
          } else if (row.annotation === 'int') {
            row.default = parseInt(row.default)
          } else if (row.annotation.startsWith("choice")) {
            row.choice = decodeChoice(row.annotation)
            row.annotation = "choice"
          } else if (row.annotation.startsWith("select")) {
            row.choice = decodeChoice(row.annotation)
            row.annotation = "select"
          }

          this.param.push(row)
        }
      }).catch((error: AxiosError) => {
        errorPrint(error)
      })
    },
    fill_path (path: string) {
      this.file = path
      axios.get(urls.file, {params: {target: path}}).then((response: any) => {
        let files: Files = {Dirs: [], Files: []}
        for (let d of response.data) {
          if (d.isdir) {
            files.Dirs.push(d)
          } else {
            if (this.postfix !== null ) {
              if (d.path.match(this.postfix)) {
                files.Files.push(d)
              }
            } else {
              files.Files.push(d)
            }
          }
        }
        this.files = files
        // this.files = response.data
      }).catch((error: any) => {
        errorPrint(error)
      })
    },
  },
  mounted() {
    this.loadParams();
    if (this.$props.func !== 'plot') {
      this.fill_path(this.file)
    }
  },
  watch: {
    func: {
      deep: true,
      handler: function() {
        this.loadParams()
      }
    }
  },
})

</script>


<style scoped>
.infinite-list {
  max-height: 300px;
  padding: 0;
  margin: 0;
  list-style: none;
}
.infinite-list .infinite-list-item {
  display: flex;
  justify-content: left;
  margin: 5px;
  color: var(--el-color-primary);
}
.infinite-list .infinite-list-item + .list-item {
  margin-top: 10px;
}
</style>
