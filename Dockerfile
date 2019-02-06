FROM debian:latest as build
RUN apt-get update -q && \
	apt-get install -y -q \
	automake \
	build-essential \
	git \
	bash-completion
RUN git clone --depth=1 https://github.com/wadedismukes/treeducken.git /treeducken
RUN cd treeducken/src/ && make install

FROM debian:latest
ENV PATH="/treeducken/:${PATH}"
COPY --from=build /treeducken /treeducken/

ENTRYPOINT ["treeducken/treeducken"]

