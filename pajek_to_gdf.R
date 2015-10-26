library(igraph);
library(NetIndices);

isInRange<-function(a,lower,higher){
	answer<-FALSE;
	if(a >= lower){
		if(a < higher){
			answer<-TRUE;
		}	
	}
	return(answer);
}

CheckColor<-function(swtl){
	auxColor<-"'100,0,0'";#super-brown
	if(isInRange(swtl,1,1.1)){
		auxColor<-"'0,255,0'";#green
	}
	if(isInRange(swtl,1.1,2.0)){
		auxColor<-"'255,255,0'";#yellow
	}
	if(isInRange(swtl,2.0,3.0)){
		auxColor<-"'255,128,0'";#orange
	}
	if(isInRange(swtl,3.0,4.0)){
		auxColor<-"'255,0,0'";#red
	}
	if(isInRange(swtl,4.0,5.0)){
		auxColor<-"'163,82,0'";#brown
	}

	return(auxColor);
}

ShortTroph<-function(net){
	matpaths<-shortest.paths(net,mode="out")
	basals<-which(degree(net,mode="in")==0)
	notbasals<-(1:length(V(net)))[-basals]
	auxSTL<-array(0,length(matpaths[1,]))
	for(i in 1:length(notbasals)){auxSTL[notbasals[i]]<-min(matpaths[basals,notbasals[i]]);}
	auxSTL<-auxSTL+1;
	return(auxSTL);
}

Generality<-function(net,L,S){
	auxGen<-array(0,S);
	for(i in 1:S){
		auxGen[i]<-degree(net,mode="in")[i]/(L/S);
	}
	return(auxGen);
}

Vulnerability<-function(net,L,S){
	auxVul<-array(0,S);
	for(i in 1:S){
		auxVul[i]<-degree(net,mode="out")[i]/(L/S);
	}
	return(auxVul);
}

OmnivoryLevel<-function(net,SWTL,S){
	matpaths<-shortest.paths(net,mode="out")
	auxOmn<-array(0,S);
	for(i in 1:S){
		preys<-which(matpaths[,i]==1);
		if(length(preys) > 1){
			auxOmn[i]<-sd(SWTL[preys]);
		}
	}
	return(auxOmn);
}

get_SWTL_omnpreyomn<-function(net,whoisomn,SWTL,Omn,SWTL_omn,web){
	auxout_a<-list();
	auxout_b<-list();
	matpaths<-shortest.paths(net,mode="out");
	for(i in 1:length(whoisomn)){
		preys_omn<-which(matpaths[,whoisomn[i]]==1);
		for(k in 1:length(preys_omn)){
			auxout_a<-c(auxout_a,SWTL_omn[i]);
			auxout_b<-c(auxout_b,SWTL[preys_omn[k]])	
			cat(paste(as.numeric(auxout_a[k]),as.numeric(auxout_b[k]),sep=","),sep="\n",file=paste("./output_TL_omn_TL_preyomn_",web,".dat",sep=""),append=TRUE);
		}
	}
	output<-cbind(as.numeric(auxout_a),as.numeric(auxout_b));
	return(output);
}

PrintPajek<-function(net,S,L,SWTL,filename){
#	for(i in 1:S){cat(paste("my_ID","my_NUMBER","my_LABEL","my_TL",sep=";"),file=filename,sep="\n");}
	cat(paste("*Vertices ",S,sep=""),file=filename,sep="\n");
	maxswtl<-max(SWTL);
	for(i in 1:S){
	
		swtl<-SWTL[i];
		swtl<-0.49*swtl/maxswtl;
		angle<-runif(1,min=0,max=2*pi);
		x<-0.5+swtl*cos(angle);			
		y<-0.5+swtl*sin(angle);			
		cat(paste(i,V(net)$id[i],x,y,sep=" "),file=filename,sep="\n",append=TRUE);
	}
	cat(paste("*Arcs",sep=""),file=filename,sep="\n",append=TRUE);
	sources<-get.edges(net,E(net))[,1];
	targets<-get.edges(net,E(net))[,2];
	for(i in 1:L){cat(paste(sources[i]+1,targets[i]+1,sep=" "),file=filename,sep="\n",append=TRUE);}
}

PrintGDF<-function(net,S,L,SWTL,filename){#prints the position and the size of the nodes
	cat("nodedef>name VARCHAR,label VARCHAR,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR",file=filename,sep="\n");
	maxswtl<-max(SWTL);
	for(i in 1:S){
		swtl<-SWTL[i];
		col<-CheckColor(swtl);
		swtl<-100*swtl/maxswtl;#SWTL [0:100]
		size<-1+swtl/25;#size of the point in the range [1:6]
		angle<-runif(1,min=0,max=2*pi);
		x<-0+(swtl*0.8)*cos(angle);#circle with points in range [-80:80]			
		y<-0+(swtl*0.8)*sin(angle);			
		cat(paste(i,V(net)$id[i],size,size,x,y,col,sep=","),file=filename,sep="\n",append=TRUE);
	}
	cat("edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,width DOUBLE,color VARCHAR",sep="\n",file=filename,append=TRUE);
	sources<-get.edges(net,E(net))[,1];
	targets<-get.edges(net,E(net))[,2];
	for(i in 1:L){cat(paste(sources[i]+1,targets[i]+1,"true",0.1,"'255,254,254'",sep=","),file=filename,sep="\n",append=TRUE);}
}

PrintGDF2<-function(net,S,L,SWTL,filename){#DOES NOT print the position and the size of the nodes
	cat("nodedef>name VARCHAR,label VARCHAR,color VARCHAR",file=filename,sep="\n");
	maxswtl<-max(SWTL);
	for(i in 1:S){
		swtl<-SWTL[i];
		col<-CheckColor(swtl);
		swtl<-100*swtl/maxswtl;#SWTL [0:100]
		size<-1+swtl/25;#size of the point in the range [1:6]
		angle<-runif(1,min=0,max=2*pi);
		x<-0+(swtl*0.8)*cos(angle);#circle with points in range [-80:80]			
		y<-0+(swtl*0.8)*sin(angle);			
		cat(paste(i,V(net)$id[i],col,sep=","),file=filename,sep="\n",append=TRUE);
	}
	cat("edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,width DOUBLE,color VARCHAR",sep="\n",file=filename,append=TRUE);
	sources<-get.edges(net,E(net))[,1];
	targets<-get.edges(net,E(net))[,2];
	for(i in 1:L){cat(paste(sources[i]+1,targets[i]+1,"true",0.1,"'255,254,254'",sep=","),file=filename,sep="\n",append=TRUE);}
}

PrintParameters<-function(net,S,L,SWTL,filename){
	cat(paste("my_ID","my_LABEL","my_TL",sep=";"),file=filename,sep="\n");
	for(i in 1:S){
		swtl<-SWTL[i];
		cat(paste(i,V(net)$id[i],swtl,sep=";"),file=filename,sep="\n",append=TRUE);
	}
}

pajek_to_gdf<-function(web){#web is the name of the file with the network in format pajek
	net<-read.graph(paste(web,sep=""),format="pajek");
	mat<-as.matrix(get.adjacency(net));
	S<-length(V(net));
	L<-length(E(net));
	C<-L/(S^2);
	basals<-which(degree(net,mode="in")==0)
	SWTL<-array(0,S);
	STL<-array(0,S);
	PATL<-array(0,S);
	PATL<-TrophInd(mat)$TL;
	STL<-ShortTroph(net);
	SWTL<-((PATL + STL)/2);
	SWTL[basals]<-1.000000;
	Omn<-OmnivoryLevel(net,SWTL,S);
	PercOmn<-100*sum(Omn>0)/S;
	whoisomn<-which(Omn > 0);
	SWTL_omn<-SWTL[whoisomn];
	SWTL_omn_preyomn<-get_SWTL_omnpreyomn(net,whoisomn,SWTL,Omn,SWTL_omn,web)

	PrintGDF(net,S,L,SWTL,paste(web,".gdf",sep=""))
	PrintGDF2(net,S,L,SWTL,paste(web,"_2.gdf",sep=""))
}

filename<-"PierreNetwork.net";
pajek_to_gdf(filename);
