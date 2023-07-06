import {AxiosError} from "axios";
import {ElNotification} from "element-plus";
import {VNode} from 'vue'

export interface Notification {
  title: string,
  message: string | VNode,
  type: string
}

export function errorPrint(error: AxiosError|Notification) {
  if (error instanceof AxiosError) {
    ElNotification({
      type: 'error',
      title: `Error Status: ${error.name}`,
      message: error.message,
      showClose: true
    })
  } else {
    ElNotification({
      type: error.type,
      title: error.title,
      message: error.message,
      showClose: true
    })
  }
}
