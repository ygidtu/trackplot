<script setup lang="ts">
import {Message, Document, Folder} from '@element-plus/icons-vue'
</script>

<template>
  <div class="params">
    <el-form ref="form" label-width="160px">
      <div v-for="(p, index) in param" :key="index">
        <el-row :gutter="20">
          <el-col :span="20">
            <el-form-item :label="p.key" v-if="!String(p.default).includes('inspect._empty')">
              <el-input v-if="p.annotation === 'str' || p.annotation === 'Optional[str]'" v-model="p.default"/>
              <el-input-number v-else-if="p.annotation === 'int'" v-model="p.default"/>
              <el-input-number v-else-if="p.annotation === 'float'" :precision='2' v-model="p.default"/>

              <el-radio-group v-model="p.default" :default="p.default" v-else-if="p.annotation === 'choice'">
                <el-radio v-for="i in p.choice" type="primary" :key="i" :label="i">{{i}}</el-radio>
              </el-radio-group>

              <el-switch v-else-if="p.annotation === 'bool'" v-model="p.default" active-text="True" inactive-text="False"/>
              <el-input v-else-if="!p.default.includes('inspect._empty')" v-model="p.default"></el-input>
            </el-form-item>
          </el-col>
          <el-col :span="4">
            <el-popover v-if="p.note !== ''"
                  placement="top-start"
                  title="Description"
                  :width="400"
                  trigger="hover"
                  :content="p.note"
                >
              <template #reference>
                <el-button :icon="Message" circle />
              </template>
            </el-popover>
          </el-col>
        </el-row>
      </div>

      <div v-if="img === null">
        <el-form-item label="File path">
          <el-input v-model="file" @input="fill_path(file)" placeholder="please select or input the file path" />
        </el-form-item>

        <el-divider />

        <ul class="infinite-list" style="overflow: auto">
          <li v-for="i in files.Dirs" :key="i.path" class="infinite-list-item" @click="fill_path(i.path)">
            <Folder />{{ i.path }}
          </li>
          <li v-for="i in files.Files" :key="i.path" class="infinite-list-item" @click="fill_path(i.path)">
            <Document />{{ i.path }}
          </li>
        </ul>

      </div>

      <el-button type="primary" @click="submit" v-if="$props.func !== 'plot'">Confirm</el-button>
      <el-button-group v-else>
        <el-button type="primary" icon="el-icon-view" @click="submit">Preview</el-button>
        <el-button type="primary" @click="save">Save pdf<i class="el-icon-download"></i></el-button>
      </el-button-group>
    </el-form>

    <div v-if="img !== null">
      <el-divider/>
      <el-row>
        <el-col :span="20" :offset="2">
          <el-image :src="img"></el-image>
        </el-col>
      </el-row>
    </div>
  </div>
</template>


<script lang="ts">
import {defineComponent} from 'vue'
import { saveAs } from 'file-saver'
import {AxiosResponse, AxiosError} from 'axios'
import urls from '../url'

import {errorPrint, Notification} from "../error.ts";


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
  img: string | null,
  param: Param[],
  file: string,
  files: Files
}

const decodeChoice = (choices: String) => {
  let res = Array<String>();
  let choiceMatch = choices.match(/choice\[(.*)\]/);
  if (choiceMatch === null) {
    return res
  }
  
  let choice = choiceMatch[1].toString().replaceAll("'", "").replaceAll('"', '').split(',');

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
      img: null, param: [],
      file: "",
      files: {
        Dirs: [],
        Files: []
      }
    }
    return d
  },
  emits: ["select-data"],
  methods: {
    save() {
      if (this.$props.func !== "plot") {
        let msg: Notification = {
          type: 'error',
          title: `Please choose file first`,
          message: ""
        }
        errorPrint(msg)
        return
      }

      this.axios.post(
          `${urls.plot}?pid=${this.$cookie.getCookie("plot")}&func=save`,
          {
            path: this.file,
            param: this.param,
          },
          {responseType: 'blob'}
      ).then((response: AxiosResponse) => {
        let headers = response.headers
        let filename = headers["content-disposition"].split("filename=")[1]
        saveAs(response.data, filename)
      }).catch((error: AxiosError) => {
        errorPrint(error)
      })
    },
    submit() {
      this.$emit("select-data", {path: this.file, param: this.param})
    },
    loadParams() {
      // 👇️ const data: GetUsersResponse
      this.axios.get(urls.params, {params: {target: this.$props.func}}
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
          }

          this.param.push(row)
        }
      }).catch((error: AxiosError) => {
        errorPrint(error)
      })
    },
    fill_path (path: string) {
      this.file = path
      this.axios.get(urls.file, {params: {target: path}}).then((response: any) => {
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
    if (this.img === null) {
      this.fill_path(this.file)
    }

  },
})

</script>


<style scoped>
.infinite-list {
  height: 300px;
  padding: 0;
  margin: 0;
  list-style: none;
}
.infinite-list .infinite-list-item {
  display: flex;
  justify-content: left;
  height: 30px;
  margin: 10px;
  color: var(--el-color-primary);
}
.infinite-list .infinite-list-item + .list-item {
  margin-top: 10px;
}
</style>