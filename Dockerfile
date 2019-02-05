FROM alpine:latest
COPY . /treeducken
RUN apt-get update && apt-get install -y \
	automake \
	curl \
	build-essential \
	git \
	bash-completion

RUN cd treeducken/src/ && make install
ENV PATH="/treeducken/:${PATH}"
RUN echo $PATH


