#include	<algorithm>
#include	<iostream>
#include	<unistd.h>
#include	<fstream>
#include	<cstring>
#include	<vector>
#include	"wyhash.h"
#include	<zlib.h>
using	namespace	std;

class	Capture{
private:
	vector<bool>	bloom;
	vector<string>	fasta;
	uint64_t	rotl(const	uint64_t	x,	int	k){	return (x << k) | (x >> (64 - k));  }
	inline	unsigned	hits(char	*S);
public:
	uint64_t	hash_len,	hash_hit;
	bool	load_fasta(const	char	*F);
	void	fasta2keys(void);
	void	capture(const	char	*O,	const	char	*R1,	const	char	*R2);
};

unsigned	Capture::hits(char	*S){
	unsigned	len=strlen(S)-1,	h=0;
	for(unsigned	j=0;	j+hash_len<=len;	j++){
		uint64_t	k=wyhash(S+j,	hash_len,	0);
		if(bloom[k&0xffffffff]&&bloom[k>>32]&&bloom[rotl(k,16)&0xffffffff])	h++;
	}	
	return	h;
}

bool	Capture::load_fasta(const	char	*F){
	ifstream	fi(F);
	if(!fi){	cerr<<"fail to open "<<F<<'\n';	return	false;	}
	string	line,	data;
	while(getline(fi,	line))	if(line.size()){
		if(line[0]=='>'){	if(data.size())	fasta.push_back(data);	data.clear();	}
		else	data+=line;
	}
	if(data.size())	fasta.push_back(data);
	fi.close();
	cerr<<"targets\t"<<fasta.size()<<'\n';
	return	true;
}

void	Capture::fasta2keys(void){
	char	fwd[256]={};	fwd['A']=fwd['a']='A';	fwd['C']=fwd['c']='C';	fwd['G']=fwd['g']='G';	fwd['T']=fwd['t']='T';	fwd['N']=fwd['n']='N';
	char	rev[256]={};	rev['A']='T';	rev['C']='G';	rev['G']='C';	rev['T']='A';	rev['N']='N';
	bloom.resize(1ULL<<32);	string	ns;	ns.assign(hash_len,	'N');	uint64_t	hashn=wyhash(ns.c_str(),	hash_len,	0);
	for(size_t	i=0;	i<fasta.size();	i++){
		string	&s=fasta[i];
		for(size_t	j=0;	j<s.size();	j++)	s[j]=fwd[(unsigned char)s[j]];
		for(size_t	j=0;	j+hash_len<=s.size();	j++){
			uint64_t	k=wyhash(s.c_str()+j,	hash_len,	0);	if(k==hashn)	continue;
			bloom[k&0xffffffff]=1;	bloom[k>>32]=1;	bloom[rotl(k,16)&0xffffffff]=1;
		}	
		reverse(s.begin(),	s.end());
		for(size_t	j=0;	j<s.size();	j++)	s[j]=rev[(unsigned char)s[j]];
		for(size_t	j=0;	j+hash_len<=s.size();	j++){
			uint64_t	k=wyhash(s.c_str()+j,	hash_len,	0);	if(k==hashn)	continue;
			bloom[k&0xffffffff]=1;	bloom[k>>32]=1;	bloom[rotl(k,16)&0xffffffff]=1;
		}
	}
	vector<string>().swap(fasta);
}

void	Capture::capture(const	char	*O,	const	char	*R1,	const	char	*R2){
	string	fn=O;	fn+="_R1.fq.gz";	gzFile	out1=gzopen(fn.c_str(),	"wt");
	fn=O;	fn+="_R2.fq.gz";	gzFile	out2=gzopen(fn.c_str(),	"wt");
	gzFile	in1=gzopen(R1,	"rt"),	in2=gzopen(R2,	"rt");
	if(in1==Z_NULL||in2==Z_NULL)	return;
	gzbuffer(in1,	1024*1024);	gzbuffer(in2,	1024*1024);
	char	nam1[256],	seq1[1024],	add1[256],	qua1[1024],	nam2[256],	seq2[1024],	add2[256],	qua2[1024];
	uint64_t	total=0,	target=0;
	while(!gzeof(in1)&&!gzeof(in2)){
		if(gzgets(in1,	nam1,	256)==NULL)	break;
		if(gzgets(in1,	seq1,	1024)==NULL)	break;
		if(gzgets(in1,	add1,	256)==NULL)	break;
		if(gzgets(in1,	qua1,	1024)==NULL)	break;	
		if(gzgets(in2,	nam2,	256)==NULL)	break;
		if(gzgets(in2,	seq2,	1024)==NULL)	break;
		if(gzgets(in2,	add2,	256)==NULL)	break;
		if(gzgets(in2,	qua2,	1024)==NULL)	break;	
		unsigned	h1=hits(seq1),	h2=hits(seq2);
		if(++total%1000000==0)	cerr<<'=';
		if(h1>=hash_hit||h2>=hash_hit){
			target++;
			gzprintf(out1,	"%s%s%s%s",	nam1,	seq1,	add1,	qua1);
			gzprintf(out2,	"%s%s%s%s",	nam2,	seq2,	add2,	qua2);
		}
	}	
	cerr<<'\n';
	gzclose(out1);	gzclose(out2);	gzclose(in1);	gzclose(in2);
	cerr<<"capture\t"<<target<<'/'<<total<<'\t'<<100.0*target/total<<"%\n";
}

void	document(void) {
	cerr<<"Usage:	iscap [options] bed.fasta output read1 read2\n";
	cerr<<"\t-l:	kmer length	default=16\n";
	cerr<<"\t-k:	#key hits	default=16\n";
	exit(0);
}

int	main(int	ac,	char	**av){
	size_t	t0=time(NULL);
	Capture	c;	c.hash_len=16;	c.hash_hit=16;
	int	opt;
	while((opt=getopt(ac,	av,	"l:k:"))>=0) {
		switch(opt) {
		case	'l':	c.hash_len=atoi(optarg);	break;
		case	'k':	c.hash_hit=atoi(optarg);	break;
		default:	document();
		}
	}
	if(ac<optind+4)	document();
	cerr<<"kmer length\t"<<c.hash_len<<'\n';
	cerr<<"#key hits\t"<<c.hash_hit<<'\n';	
	if(!c.load_fasta(av[optind]))	return	0;
	c.fasta2keys();
	c.capture(av[optind+1],	av[optind+2],	av[optind+3]);
	cerr<<"time\t"<<time(NULL)-t0<<" s\n";
	return	0;
}
