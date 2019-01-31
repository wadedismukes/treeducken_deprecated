FROM debian:latest
COPY . /treeducken
RUN apt-get update && apt-get install -y \
	automake \
	curl \
	build-essential \
	git \
	curl \
	bash-completion

RUN cd treeducken/src/ && make install
ENV PATH="/treeducken/:${PATH}"


