# ==================== 第一阶段：构建依赖 ====================
FROM python:3.12-slim AS builder

# 安装编译依赖（部分包需要编译）
# 默认使用官方源 deb.debian.org（国外用户开箱即用）；国内用户可 --build-arg MIRROR=mirrors.aliyun.com 覆盖
ARG MIRROR=
RUN if [ -n "${MIRROR}" ]; then \
        find /etc/apt \( -name "*.list" -o -name "*.sources" \) -print0 2>/dev/null | xargs -0 -r sed -i -e "s|deb.debian.org/debian|${MIRROR}/debian|g" -e "s|security.debian.org/debian-security|${MIRROR}/debian-security|g" 2>/dev/null || true; \
    fi && \
    apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# 配置 pip 使用本地源（可选）
ARG PIP_INDEX_URL=https://pypi.org/simple
ARG PIP_TRUSTED_HOST=pypi.org

WORKDIR /build

# 只复制依赖文件，利用 Docker 缓存
COPY pyproject.toml requirements.txt ./
COPY ./trackplot ./trackplot

# 构建 wheel 包
# 创建 .local 目录并配置 pip
RUN mkdir -p /root/.local/bin && \
    pip config set global.index-url $PIP_INDEX_URL && \
    if [ "$PIP_TRUSTED_HOST" != "pypi.org" ]; then \
        pip config set global.trusted-host $PIP_TRUSTED_HOST; \
    fi && \
    # 升级 pip 和安装工具
    pip install --no-cache-dir --upgrade pip setuptools wheel --default-timeout=100 --retries 5 && \
    # 安装项目（使用 --prefix 确保文件在已知位置）
    pip install --no-cache-dir --prefix=/root/.local -e . --default-timeout=100 --retries 5 || \
    pip install --no-cache-dir --target=/root/.local/lib/python3.11/site-packages -e . --default-timeout=100 --retries 5

# 验证安装
RUN ls -la /root/.local/bin/ || echo "No binaries installed"


# ==================== 第二阶段：运行环境 ====================
FROM python:3.12-slim

ARG PIP_INDEX_URL=https://pypi.org/simple
ENV PIP_INDEX_URL=$PIP_INDEX_URL

# 安装运行时依赖（只安装必要的系统库）
# 默认使用官方源 deb.debian.org（国外用户开箱即用）；国内用户可 --build-arg MIRROR=mirrors.aliyun.com 覆盖
ARG MIRROR=
RUN if [ -n "${MIRROR}" ]; then \
        find /etc/apt \( -name "*.list" -o -name "*.sources" \) -print0 2>/dev/null | xargs -0 -r sed -i -e "s|deb.debian.org/debian|${MIRROR}/debian|g" -e "s|security.debian.org/debian-security|${MIRROR}/debian-security|g" 2>/dev/null || true; \
    fi && \
    apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 从构建阶段复制已安装的包
COPY --from=builder /root/.local /root/.local

ENV PATH=/root/.local/bin:$PATH
ENV ROOT_DIR=/opt/trackplot

WORKDIR $ROOT_DIR
COPY ./ $ROOT_DIR

# 安装本地包（将 trackplot 安装到系统 Python，源文件在 /opt/trackplot）
# 加 --retries/--default-timeout 应对慢网络；PIP_INDEX_URL 已在国内镜像下生效
RUN pip install --no-cache-dir --default-timeout=100 --retries 5 -e .

ENTRYPOINT ["python", "/opt/trackplot/main.py"]