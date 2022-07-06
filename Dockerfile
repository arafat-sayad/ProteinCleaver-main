FROM gkoulouras/proteincleaver-packages:1.0.0

ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    sudo \
    git 
RUN R -q -e "devtools::install_github('sgibb/cleaver')"

   
RUN mkdir /root/app


COPY ./src /root/app/src


COPY Rprofile.site /usr/lib/R/etc/


EXPOSE 3838

# CMD Rscript /root/app/src/app.R

CMD ["R", "-e", "shiny::runApp('/root/app/src')"]

# sudo docker build --tag=gkoulouras/proteincleaver-application:1.0.0 -f Dockerfile .
# sudo docker run  --name bi --publish 3838:3838 gkoulouras/proteincleaver-application:1.0.0 
