<template>

  <el-descriptions>
    <el-descriptions-item>
      <el-text type="info">If there is any unexpected error, please send </el-text>
      <el-link type="primary" :href="`${urls.log}?pid=${pid}&download=true`">this log file</el-link>
      <el-text type="info"> to developers at </el-text>
      <el-link href="https://github.com/ygidtu/trackplot/issues" type="primary">Github</el-link>
      <el-text type="info"> for debugging</el-text>
    </el-descriptions-item>
  </el-descriptions>
  <el-divider />
    <el-timeline>
      <el-timeline-item
        v-for="activity in logs"
        :type="logLevel(activity.level)"
        :key="activity.time"
        :size="'large'"
        :center="false"
        placement="top"
        hide-timestamp
      >
        <el-descriptions>
          <el-descriptions-item :label="activity.time">
            <el-text tag="b">{{ activity.source }}</el-text>
            <el-divider direction="vertical" />
            <el-text :type="logLevel(activity.level)">{{ activity.message }}</el-text>
          </el-descriptions-item>
        </el-descriptions>
      </el-timeline-item>
    </el-timeline>
</template>

<script lang="ts" setup>
import urls from "../url";
</script>

<script lang="ts">
import {AxiosResponse, AxiosError} from "axios";
import {errorPrint} from "../error";
interface Log {
  time: string,
  level: string,
  source: string,
  message: string
}
export default {
  name: 'LogComp',
  props: {
    pid: {required: false, type: String, default: "test"},
    load: {required: true, type: Boolean, default: true}
  },
  data() {
    let logs: Array<Log> = []
    return {
      timer: "",
      logs: logs,
      count: 1
    }
  },
  methods: {
    loadParams() {
      // 👇️ const data: GetUsersResponse
      this.axios.get(urls.log, { params: {pid: this.$props.pid}}
      ).then((response: AxiosResponse) => {
        this.logs = response.data
      }).catch((error: AxiosError) => {
        if (this.count < 2) {
          errorPrint(error)
          this.count = 2
        }
      })
    },
    logLevel (level: string) {
      if (level === "INFO") {
        return("primary")
      }
      if (level === "WARN") {
        return("warning")
      }
      if (level === "ERROR") {
        return("danger")
      }
      return("info")
    }
  },
  mounted() {
    this.timer = setInterval(this.loadParams, 1000);
  },
  watch: {
    load: function() {
      if (this.load) {
        this.timer = setInterval(this.loadParams, 1000);
      } else {
        clearInterval(this.timer);
      }
    }
  },
  beforeUnmount() {
    clearInterval(this.timer);
  }
}
</script>

<style scoped>

</style>