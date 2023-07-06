<template>

  <el-descriptions>
    <el-descriptions-item :label="'If there is any unexpected error, please send this log file to developers for debugging'">
      <el-link type="primary" :href="`${urls.log}?pid=${pid}&download=true`">Download</el-link>
    </el-descriptions-item>
  </el-descriptions>
  <el-divider />
  <el-scrollbar :height="height" always>
    <el-timeline>
      <el-timeline-item
        v-for="activity in logs"
        :type="this.logLevel(activity.level)"
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
            <el-text :type="this.logLevel(activity.level)">{{ activity.message }}</el-text>
          </el-descriptions-item>
        </el-descriptions>
      </el-timeline-item>
    </el-timeline>
  </el-scrollbar>
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
      height: {default: "500px"}
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
          errorPrint(error)
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
      // this.timer = setInterval(this.loadParams, 2000);
    },
    beforeUnmount() {
      clearInterval(this.timer);
    }
  }
</script>

<style scoped>

</style>