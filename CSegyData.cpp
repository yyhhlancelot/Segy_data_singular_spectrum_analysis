#include "stdafx.h"
#include <math.h>
#include  <fcntl.h>
#include <iostream>
#include "CSegyData.h"
using namespace std;

CSegyData::CSegyData()
{
	m_segyFileName = NULL;
	m_buffer = NULL ;
	m_bSegyFileOpened = false;
	m_fileSegy = 0;
}

CSegyData::CSegyData(char* segyFileName,int segyFileType)
{
	OpenSegyFile(segyFileName,segyFileType);
}

CSegyData::~CSegyData()
{
	CloseSegyFile( );
}


/************************************************************
Function:       ReadOneTrace(__int64 nTraceIndex, TraceHead *pTraceHead, float *pTraceData)
Description:    根据输入道序号，读取一个地震道
Calls:          无
Inputs:         输入道序号，从0开始计数( 范围： 0 <--> GetSegyTotalTraceNumber() - 1 )、道头指针、地震数据数组指针（长度为GetSegySampleNumber()）地震数据数组指针（长度为GetSegySampleNumber()）
Output:         无
Return:         无 
Others:         返回值：是否读取成功，通常，如果输入道序号范围正确，总是会返回true
*************************************************************/
bool CSegyData::ReadOneTrace(__int64 nTraceIndex, TraceHead *pTraceHead, float *pTraceData)
{
	if (!m_bSegyFileOpened)
	{
		return false;
	}

	if (nTraceIndex< 0 || m_fileSegy >= GetSegyTotalTraceNumber())
	{
		return false;
	}

	__int64 offset = nTraceIndex * m_traceBytes + 3200 + 400 ;
	_lseeki64(m_fileSegy, offset, SEEK_SET);

	if( read( m_fileSegy, m_buffer , m_traceBytes ) != m_traceBytes )
	{
		printf("File %s read %d bytes error ! \n", m_segyFileName,m_traceBytes );
		return false;
	}

	unsigned long  * pInt32 =(unsigned long  * )&m_buffer[240];

	if( m_segyFileType == 0)  //0 --- tape file  交换字节顺序
	{
		TraceHead *pTraceHead = (TraceHead *) m_buffer;
		swapSegyTraceHeadBytes( *pTraceHead );//240 bytes
		SwapInt4(&m_buffer[240], GetSegySampleNumber());
	}
	else  // 1----disk file 不用交换字节顺序
	{
	}

	if(m_lineHead.dataform == 1) //解编为标准的二进制浮点格式
	{
		int k;
		for(k = 0 ; k < m_lineHead.samplenum ; k++ )
		{
			if(pInt32[k] != 0 )
			{
				pInt32[k] = SegYToFloat(  pInt32[k] );
			}
		}
	}

	memcpy(pTraceHead,  m_buffer,      240);
	memcpy(pTraceData, &m_buffer[240], m_lineHead.samplenum * 4);
	return true;
}

/************************************************************
Function:       ReadOneTrace(__int64 nTraceIndex, float *pTraceData)
Description:    根据输入道序号，读取一个地震道
Calls:          无
Inputs:         输入道序号，从0开始计数( 范围： 0 <--> GetSegyTotalTraceNumber() - 1 )、地震数据数组指针（长度为GetSegySampleNumber()）地震数据数组指针（长度为GetSegySampleNumber()）
Output:         无
Return:         无 
Others:         返回值：是否读取成功，通常，如果输入道序号范围正确，总是会返回true
*************************************************************/
bool CSegyData::ReadOneTrace(__int64 nTraceIndex, float *pTraceData)
{
	if (!m_bSegyFileOpened)
	{
		return false;
	}

	if (nTraceIndex< 0 || m_fileSegy >= GetSegyTotalTraceNumber())
	{
		return false;
	}

	__int64 offset = nTraceIndex * m_traceBytes + 3200 + 400 ;
	_lseeki64(m_fileSegy, offset, SEEK_SET);

	if( read( m_fileSegy, m_buffer , m_traceBytes ) != m_traceBytes )
	{
		printf("File %s read %d bytes error ! \n", m_segyFileName,m_traceBytes );
		return false;
	}

	unsigned long  * pInt32 =(unsigned long  * )&m_buffer[240];
	if( m_segyFileType == 0 )  //0 --- tape file  交换字节顺序
	{
		TraceHead *pTraceHead = (TraceHead *) m_buffer;
		SwapInt4(&m_buffer[240], GetSegySampleNumber());
	}
	else  // 1----disk file 不用交换字节顺序
	{
	}

	if(m_lineHead.dataform==1) //解编为标准的二进制浮点格式
	{
		int k;
		for(k=0 ; k < m_lineHead.samplenum ; k++ )
		{
			if(pInt32[k] != 0 )
			{
				pInt32[k] = SegYToFloat(  pInt32[k] ) ;
			}
		}
	}

	memcpy(pTraceData, &m_buffer[240], m_lineHead.samplenum * 4);

	return true;
}

/************************************************************
Function:       ReadOneTraceHeader(__int64 nTraceIndex, TraceHead *pTraceHead)
Description:    根据输入道序号，读取一个地震道的道头
Calls:          无
Inputs:         输入道序号，从0开始计数( 范围： 0 <--> GetSegyTotalTraceNumber() - 1 )、道头指针
Return:         无 
Others:         返回值：是否读取成功，通常，如果输入道序号范围正确，总是会返回true
*************************************************************/
bool CSegyData::ReadOneTraceHeader(__int64 nTraceIndex, TraceHead *pTraceHead)
{
	if (!m_bSegyFileOpened)
	{
		return false;
	}

	if (nTraceIndex< 0 || m_fileSegy >= GetSegyTotalTraceNumber())
	{
		return false;
	}

	__int64 offset = nTraceIndex * m_traceBytes + 3200 + 400 ;
	_lseeki64(m_fileSegy, offset, SEEK_SET);

	if( read( m_fileSegy, m_buffer , 240 ) != 240 )
	{
		printf("File %s read traces %d, %d bytes error ! \n", m_segyFileName,nTraceIndex,m_traceBytes );
		return false;
	}

	if( m_segyFileType == 0 )  //0 --- tape file  交换字节顺序
	{
		TraceHead *pTraceHead = (TraceHead *) m_buffer;
		swapSegyTraceHeadBytes( *pTraceHead ); //240 bytes
	}
	else  // 1----disk file 不用交换字节顺序
	{
	}

	memcpy(pTraceHead, m_buffer, 240);
	return true;
}

/************************************************************
Function:       SplitSegyFile(const char *newSegyFileName ,const __int64 m_jumpSamp,const __int64 m_jumpTrace,
                const  short m_showSamp,const  int 		m_showTrace	)
Description:    把原始segy文件提取一部分放到新的segy文件里面
Calls:          无
Inputs:         新segy文件的存储地址、m_jumpSamp: 每一道跳过的样点数、m_showSamp: 每一道显示的样点数、m_jumpTrace: 从文件头要跳过的道数、m_showTrace:要截取的道数
Return:         无 
Others:         无
*************************************************************/
bool CSegyData::SplitSegyFile(const char *newSegyFileName ,
							 const __int64 		m_jumpSamp,
							 const __int64 		m_jumpTrace,
							 const  short 		m_showSamp,
							 const  int 		m_showTrace	)
{
	if( m_showSamp <=0 || m_showTrace <= 0 || 
		m_jumpSamp < 0 || m_jumpTrace< 0 || 
		(m_showSamp > (GetSegySampleNumber() - m_jumpSamp )) ||
		(m_showTrace> (GetSegyTotalTraceNumber()-m_jumpTrace))
		)
	{
        printf("SegyToBinary()参数中含有无效的值, 请改正 !\n");
		return false;
	}
	int fdi = open((const char *)m_segyFileName,O_RDONLY|O_BINARY);
	if( fdi == NULL ) 
	{
		printf("File %s open error ! \n", m_segyFileName );
		return false;
	}

	int fdo = open((const char *)newSegyFileName,O_WRONLY|O_TRUNC|O_CREAT|O_BINARY,0644);
	if( fdo == NULL ) 
	{
		printf("File %s open error ! \n", newSegyFileName );
		return false;
	}
    //分割文件的文本头还是原样复制
	int i;
	char * buffer	;
	buffer = new char [3200];
	i = read(fdi,buffer ,3200);
	if(i!=3200)
	{
		printf("File %s read 3200 bytes error ! \n",m_segyFileName);
		return 1;
	}
	i= write(fdo,buffer ,3200);
	if(i!=3200)
	{
		printf("File %s write 3200 bytes error ! \n",m_segyFileName);
		return false;
	}
	//分割的二进制头，
	i= read(fdi,buffer ,400);
	if(i!=400)
	{
		printf("File %s read 400 bytes error ! \n",m_segyFileName);
		return false;
	}
	//需要改变样点数值，其他二进制成员变量根据需要来变		   
	((LineHead  *)buffer)->samplenum= m_showSamp;
	((LineHead  *)buffer)->samplenum0= m_showSamp ; 
	if( m_segyFileType == 0 )
	{
		//((LineHead  *)buffer)->samplenum= SwapInt2( m_showSamp);
		//((LineHead  *)buffer)->samplenum0= SwapInt2( m_showSamp) ; 
		((LineHead  *)buffer)->samplenum = 
			SwapInt2( ((LineHead  *)buffer)->samplenum );
		((LineHead  *)buffer)->samplenum0 = 
			SwapInt2( ((LineHead*)buffer)->samplenum0 ) ;
	}
	i= _write(fdo,buffer ,400);
	if(i!=400)
	{
		printf("File %s write 400 bytes error ! \n",m_segyFileName);
		return 1;
	}

	__int64 offset = m_jumpTrace * m_traceBytes ;
	_lseeki64(fdi,offset ,SEEK_CUR);

	int outTraceBytes,skipBytes ;	
	if( m_lineHead.dataform==1||
		m_lineHead.dataform==2||
		m_lineHead.dataform==4||
		m_lineHead.dataform==5 )
	{
		outTraceBytes = m_showSamp * 4 + 240 ;
		skipBytes = m_jumpSamp *4 + 240 ;
		//skipBytes = m_jumpSamp *4 ;
	}
	else if(m_lineHead.dataform==3)
	{
		outTraceBytes = m_showSamp * 2 + 240 ;
		skipBytes = m_jumpSamp *2 + 240 ;
		//skipBytes = m_jumpSamp *2 ;
	}

	delete buffer;
	buffer = new char [m_traceBytes];

	for(long j=1;j<=m_showTrace;j++)
	{
		if (j%1000 ==0)
		{
			printf("分割完成比率 %.3f\n", ((j-1.0)/m_showTrace) *100);
		}

		if( read( fdi, buffer , m_traceBytes ) != m_traceBytes )
		{
			printf("File %s read %d bytes error ! \n", m_segyFileName,m_traceBytes );
			break ;
		}
	
		if( m_segyFileType == 0 )  //0 --- tape file  1----disk file
		{	
		/*	((TraceHead  *)buffer)->NumberSamples = 
				SwapInt2(((TraceHead  *)buffer)->TrNumCdp );*/
			
			((TraceHead  *)buffer)->NumberSamples = SwapInt2( m_showSamp) ;

			((TraceHead *)buffer)->TraceLsNum = SwapInt4(j);
			((TraceHead *)buffer)->TraceRsNum = SwapInt4(j);
		}
		else
		{
			((TraceHead *)buffer)->TraceLsNum = j;
			((TraceHead *)buffer)->TraceRsNum = j;
		}

		//buffer 里面的值本来就是segy 格式的不需要转换
		memcpy( buffer + 240 , buffer + skipBytes , outTraceBytes -240);

		if( write( fdo, buffer , outTraceBytes ) != outTraceBytes )
		{
			printf("File %s write %d bytes error ! \n", newSegyFileName,outTraceBytes );
			break ;
		}
	}

	delete buffer;
	close(fdi);
	close(fdo);
	return 0;
}

/************************************************************
Function:       SegyToBinary(char * w_buffer,const int &m_jumpSamp,const int &m_jumpTrace,
                const int &m_showSamp,const int &m_showTrace	)
Description:    把segy文件截取并保存到缓冲区里面，w_buffer不包括道头信息
Calls:          无
Inputs:         缓冲区w_buffer、m_jumpSamp: 每一道跳过的样点数、m_showSamp: 每一道显示的样点数、m_jumpTrace: 从文件头要跳过的道数、m_showTrace:要截取的道数
Return:         无 
Others:         w_buffer要足够长
*************************************************************/
bool CSegyData::SegyToBinary(char * w_buffer,//its length must long enough
							 const int &m_jumpSamp,
							 const int &m_jumpTrace,
							 const int &m_showSamp,
							 const int &m_showTrace	)
{
	if( m_showSamp <=0 || m_showTrace <= 0 || 
		m_jumpSamp < 0 || m_jumpTrace< 0 || 
		(m_showSamp > (GetSegySampleNumber() - m_jumpSamp )) ||
		(m_showTrace> (GetSegyTotalTraceNumber()-m_jumpTrace))
		)
	{
		cout << " segy to binary in buffer failed" <<endl;
		return 1;
	}
	int outTraceBytes; //每一道要显示的样点数的字节数
	int skipBytes ;	//每一道要跳过的字节数

	if( m_lineHead.dataform==1||
		m_lineHead.dataform==2||
		m_lineHead.dataform==4||
		m_lineHead.dataform==5 )
	{
		outTraceBytes = m_showSamp * 4  ;
		skipBytes = m_jumpSamp *4 + 240 ;
	}
	else if(m_lineHead.dataform==3)
	{
		outTraceBytes = m_showSamp * 2 ;
		skipBytes = m_jumpSamp *2 + 240 ;
	}

	if ( w_buffer == NULL )
	{
		return outTraceBytes * m_showTrace + sizeof(int) * m_showTrace * 3 ; //data & CDP & Tracersnum & Tracerlnum value for per trace
	}

	int fdi = open( (const char *) m_segyFileName,O_RDONLY|O_BINARY);
	if( fdi == NULL ) 
	{
		cout <<"File open failed ! " << endl;
		return 1;
	}
    //定位到要跳过的道的后面
	__int64 offset = m_jumpTrace * m_traceBytes + 3200 + 400 ;
	_lseeki64(fdi,offset ,SEEK_SET);

	unsigned long  * pInt32 = nullptr;
	unsigned short * pInt16 = nullptr;
	long j,k;
	char * pp = w_buffer;//数据存储到二进制文件，w_buffer，
	char * buffer = new char [m_traceBytes];//存一个道的数据
	//int * cdplist=  new int [m_showTrace] ;//每一个道对应的cdp号
	//int * tracersnum=  new int [m_showTrace] ;//每一个道对应的线内道序列号
	//int * tracelsnum=  new int [m_showTrace] ;//每一个道对应的卷道序列号

	for(j=0;j<m_showTrace;j++)
	{//从新定位的位置开始读取，一次读取一个道
		if( read( fdi, buffer , m_traceBytes ) != m_traceBytes )
		{
			cout <<" Segy to binary  open file failed" <<endl;
			break ;
		}

		if( m_lineHead.dataform==1||
			m_lineHead.dataform==2||
			m_lineHead.dataform==4||
			m_lineHead.dataform==5 )
			pInt32 =(unsigned long  * )(buffer + skipBytes);
		else if(m_lineHead.dataform==3)
			pInt16 =(unsigned short * )(buffer + skipBytes);

		if( m_segyFileType == 0 )  //0 --- tape file  交换字节顺序
		{
			//cdplist[j] = SwapInt4( ((TraceHead *)buffer)->CdpNum ) ;
			//tracersnum[j] = SwapInt4( ((TraceHead *)buffer)->TraceRsNum ) ;
			//tracelsnum[j] = SwapInt4( ((TraceHead *)buffer)->TraceLsNum ) ;
			if(m_lineHead.dataform==3)
			{
				for(k=0 ; k < m_showSamp ; k++ )
					if ( pInt16[k] != 0 ) pInt16[k] = SwapInt2( pInt16[k] )  ;
			}
			else
			{//每个转换本道的每个样点
				for(k=0 ; k < m_showSamp ; k++ )
					if( pInt32[k] != 0 ) pInt32[k] = SwapInt4( pInt32[k] )  ;
			}
		}
		else  // 1----disk file 不用交换字节顺序
		{
			//cdplist[j] = ((TraceHead *)buffer)->CdpNum  ;
			//tracersnum[j] =  ((TraceHead *)buffer)->TraceRsNum  ;
			//tracelsnum[j] =  ((TraceHead *)buffer)->TraceLsNum  ;
		}

		if(m_lineHead.dataform==1) //解编为标准的二进制浮点格式
		{
			for(k=0 ; k < m_showSamp ; k++ )
			{
				if(pInt32[k] != 0 )
					pInt32[k] = SegYToFloat(  pInt32[k] ) ;
			}
		}
      //从本道样点出复制 要抽取的道数的字节数
		memcpy ( pp , buffer + skipBytes , outTraceBytes );
		pp +=  outTraceBytes ;
	}
	//上面的循环已经把要抽取的道的数据顺序存放到w_buffer里面了
	//

	delete buffer;
	close(fdi);

	return 0;
}//end SegyToBinary()2

bool CSegyData::SegyToBinary(
							const __int64 		m_jumpSamp,
							const __int64		m_jumpTrace,
							const int &		m_showSamp,
							const int &		m_showTrace,
							const char * newBinFileName)
{
	if( m_showSamp <=0 || m_showTrace <= 0 || 
		m_jumpSamp < 0 || m_jumpTrace< 0 || 
		(m_showSamp > (GetSegySampleNumber() - m_jumpSamp )) ||
		(m_showTrace> (GetSegyTotalTraceNumber()-m_jumpTrace))
		)
	{
	  printf("SegyToBinary()参数中含有无效的值, 请改正 !\n");
		return false;
	}

	int outTraceBytes,skipBytes ;	

	if( m_lineHead.dataform==1||
		m_lineHead.dataform==2||
		m_lineHead.dataform==4||
		m_lineHead.dataform==5 )
	{
		outTraceBytes = m_showSamp * 4  ;
		skipBytes = m_jumpSamp *4 + 240 ;
	}
	else if(m_lineHead.dataform==3)
	{
		outTraceBytes = m_showSamp * 2 ;
		skipBytes = m_jumpSamp *2 + 240 ;
	}

	int fdi = open((const char *)m_segyFileName,O_RDONLY|O_BINARY);
	if( fdi == NULL ) 
	{
		printf("File %s open error ! \n", m_segyFileName );
		return false;
	}

	int fdo = open((const char *)newBinFileName,O_WRONLY|O_TRUNC|O_CREAT|O_BINARY,0644);
	if( fdo == NULL ) 
	{
		printf("File %s open error ! \n", newBinFileName );
		return false;
	}

	__int64 offset = m_jumpTrace * m_traceBytes + 3200 + 400 ;
	_lseeki64(fdi,offset ,SEEK_SET);

	unsigned long  * pInt32 = nullptr;
	//long  * pInt32;
	//float * pInt32;
	unsigned short * pInt16 = nullptr;
	long j,k;
	char * buffer = new char [m_traceBytes];
	int * cdplist = new int [m_showTrace] ;
	int * tracersnum=  new int [m_showTrace] ;
	int * tracelsnum=  new int [m_showTrace] ;

	for( j=0;j < m_showTrace;++j)
	{
    
       if (j %1000 == 0)
       {
		   printf("已写入的道数 %d, 完成比率 %.3f\n", j, (j +1.0)/m_showTrace *100);
       }
          

		if( read( fdi, buffer , m_traceBytes ) != m_traceBytes )
		{
			printf("File %s read %d bytes error ! \n", m_segyFileName,m_traceBytes );
			break ;
		}

		if( m_lineHead.dataform==1||
			m_lineHead.dataform==2||
			m_lineHead.dataform==4||
			m_lineHead.dataform==5 )
			pInt32 =(unsigned long  * )(buffer + skipBytes);

		else if(m_lineHead.dataform==3)
			pInt16 =(unsigned short * )(buffer + skipBytes);

		if( m_segyFileType == 0 )  //0 --- tape file  交换字节顺序
		{
			//cdplist[j] = SwapInt4( ((TraceHead *)buffer)->cdpensemblenum ) ;
			//tracersnum[j] = SwapInt4( ((TraceHead *)buffer)->tracersnum ) ;
			//tracelsnum[j] = SwapInt4( ((TraceHead *)buffer)->tracelsnum ) ;
			if(m_lineHead.dataform==3)
			{
				for(k=0 ; k < m_showSamp ; k++ )
					if ( pInt16[k] != 0 ) pInt16[k] = SwapInt2( pInt16[k] )  ;
			}
			else
			{
				for(k=0 ; k < m_showSamp ; k++ )
					if( pInt32[k] != 0 )  pInt32[k] = SwapInt4( pInt32[k] )  ;
			}
		}
		else  // 1----disk file 不用交换字节顺序
		{
			//cdplist[j] = ((TraceHead *)buffer)->cdpensemblenum  ;
			//tracersnum[j] =  ((TraceHead *)buffer)->tracersnum  ;
			//tracelsnum[j] =  ((TraceHead *)buffer)->tracelsnum  ;
		}

		if(m_lineHead.dataform==1) //解编为标准的二进制浮点格式
		{
			for(k=0 ; k < m_showSamp ; k++ )
			{
				if(pInt32[k] != 0 )
					pInt32[k] = SegYToFloat(  pInt32[k] ) ;
			}
		}

		if( write( fdo, buffer + skipBytes , outTraceBytes ) != outTraceBytes )
		{
			printf("File %s read %d bytes error ! \n", (const char *)newBinFileName,outTraceBytes );
			break ;
		}
	}
	delete buffer;
	close(fdi);
	close(fdo);

	
	delete	cdplist ;
	delete	tracelsnum;
	delete	tracersnum;

	return 0;
}


bool CSegyData::BinaryToSegy(const char *w_BINFileName,
					         const int &w_binDataType,  //1  --- 32bit float, 2----32bit integer 3----16bit integer
							 const int &w_segyFileType, //0 --- tape file  1----disk file
							 const char   *newSegyFileName,
							 const int &w_sampNum,
							 const int &w_startCDP,
							 const int &w_traceNum,
							 const int &w_timer) //时间上抽样间隔（1,2,3……）
	{
		LineHead	lineHead; //400 bytes
		TraceHead	traceHead; //240 bytes

		memset( (char *)&lineHead,  0, 400);
		memset( (char *)&traceHead, 0, 240);

		/********************************************************************************
		Line Head Messages 
		********************************************************************************/
		lineHead.jobnum = 1; //一般没有值
		lineHead.linenum= 1 ; //线号从二进制读进来也一般没有值
		lineHead.reelnum= 1 ; //卷号一般也没有值	
		lineHead.tracenum= 0; //number of data traces per record 
		lineHead.auxtrace= 0 ;//number of auxiliary traces per record 
		lineHead.interval= w_timer * 1000 ; //for this reel of data抽样间隔
		lineHead.interval0=w_timer * 1000;	//for original field recording 
		lineHead.samplenum= w_sampNum;      //for this reel of data 抽样数
		lineHead.samplenum0= w_sampNum ;    //for original field recording

		//1 = float point(4B),2 = fixed point(4B) 3 = fixed point(2B)	4 = fixed point w/gain code(4B)          
		lineHead.dataform= w_binDataType ;  
		lineHead.cdpfold= 0 ;

		/********************************************************************************
		Trace Head Messages 
		********************************************************************************/
		/* trace sequence number within line -- */
		/* number continue to increase if additional */
		/* reels are required on same line */
		traceHead.TraceLsNum = 1;
		/* trace sequence number within reel -- */
		/* each reel stars with trace number one */
		traceHead.TraceRsNum = 1;
		/* original field record number */
		traceHead.RecordNum = 0 ;
		/* within original field record */
		traceHead.TrNum = 0;
		/*CDP ensemble number */
		traceHead.CdpNum= w_startCDP ;
		/* within the CDP ensemble -- each */
		/* ensemble stars with trace number one */
		traceHead.TrNumCdp = 0;//每个cdp下的道数
		/* 1 = seismic data, etc. */
		traceHead.TraceId = 1;
		/* 1 = production , 2 = test */
		traceHead.DataUse =  1;
		/*  in this trace */
		traceHead.NumberSamples = w_sampNum  ;
		/*  in this trace */
		traceHead.SampleInterval= w_timer * 1000 ;

		//如果是磁带seg-y文件，需要交换各个头结构的字节顺序
		if( w_segyFileType == 0 )  //0 --- tape file  1----disk file
		{
			swapSegyLineHeadBytes( lineHead); //400 bytes
			swapSegyTraceHeadBytes(traceHead ); //240 bytes
		}

		int fdi = open((const char *)w_BINFileName,O_RDONLY|O_BINARY);
		if( fdi == NULL ) 
		{
			printf("File %s open error ! \n",(const char *)w_BINFileName);
			return 1;
		}

		int fdo = open((const char *)newSegyFileName,O_WRONLY|O_TRUNC|O_CREAT|O_BINARY,0644);
		if( fdo == NULL ) 
		{
			printf("File %s open error ! \n",(const char *)newSegyFileName);
			return 1;
		}

		int i;
		char * buffer,*p;

		//write 3200 EBCDIC code 文本头用@码代替不管
		p = buffer = new char [3200];
		for(i=0;i<3200;i++)*p++='@';
		i= write(fdo,buffer ,3200);
		delete buffer;
		if(i!=3200)
		{
			printf("File %s write 3200 bytes error ! \n",(const char *)newSegyFileName);
			return 1;
		}

		i= write(fdo,(char *)&lineHead,400);
		if(i!=400)
		{
			printf("File %s write 400 bytes error ! \n",(const char *)newSegyFileName);
			return 1;
		}


		int outTraceBytes , inTraceBytes ;

		if( w_binDataType==1 ||  w_binDataType==2 )
		{
			inTraceBytes  = w_sampNum  * 4 ;
			outTraceBytes = w_sampNum  * 4 + 240 ;
		}
		else if(w_binDataType==3)
		{
			inTraceBytes  = w_sampNum  * 2 ;
			outTraceBytes = w_sampNum  * 2 + 240 ;
		}


		buffer = new char [outTraceBytes];

		unsigned long  *pInt32 = nullptr;
		unsigned short *pInt16 = nullptr;
		int k;
//按照道数来读取
		for(long j=0;j<w_traceNum;j++)
		{
			if (j %1000 == 0)
			{
				printf("完成道数 %d, 比率 %.3f\n", j, (j-1.0)/w_traceNum*100);
			}
			//如果是磁带seg-y文件，需要交换各个头结构的字节顺序
			if( w_segyFileType == 0 )  //0 --- tape file  1----disk file
			{
				//头数也要变
				traceHead.TraceLsNum = traceHead.TraceRsNum = SwapInt4 ( j+1 );
				traceHead.CdpNum= SwapInt4 ( w_startCDP + j ) ;
			}
			else
			{
				traceHead.TraceLsNum = traceHead.TraceRsNum = j+1 ;
				traceHead.CdpNum=  w_startCDP + j  ;
			}

			memcpy(buffer,(char *)&traceHead,240);

			if( read( fdi, buffer+240 , inTraceBytes ) != inTraceBytes )
			{
				printf("File %s read %d bytes error ! \n",(const char *)w_BINFileName,inTraceBytes );
				break ;
			}

			if( w_binDataType==1 || w_binDataType==2 ) 
			{
				pInt32 =(unsigned long  * )(buffer + 240);
				for(k=0 ; k < w_sampNum  ; k++ )//解编为标准的二进制浮点格式
					if(pInt32[k] != 0 )pInt32[k] = FloatToSegy(  pInt32[k] ) ;
			}
			else if(w_binDataType==3)
				pInt16 =(unsigned short * )(buffer + 240);

			if( w_segyFileType == 0 )  //0 --- tape file  交换字节顺序
			{
				if(w_binDataType==3)
				{
					for(k=0 ; k < w_sampNum ; k++ )
						if ( pInt16[k] != 0 ) pInt16[k] = SwapInt2( pInt16[k] )  ;
				}
				else
				{
					for(k=0 ; k < w_sampNum ; k++ )
						if( pInt32[k]  != 0 ) pInt32[k] = SwapInt4( pInt32[k] )  ;
				}
			}

			if( write( fdo, buffer , outTraceBytes ) != outTraceBytes )
			{
				printf("File %s read %d bytes error ! \n",(const char *)newSegyFileName,outTraceBytes );
				break ;
			}
		}

		delete buffer;
		close(fdi);
		close(fdo);
		return 0;
	}

bool CSegyData::OpenSegyFile(const char * segyFileName,const int segyFileType)
{
	//如果先前打开过一个SEGY文件，首先关闭
	if (m_bSegyFileOpened)
	{
		CloseSegyFile();
	}

	m_buffer = NULL ;
	m_bSegyFileOpened = false;

	//判断该文件是否存在
	if(_access((const char *)segyFileName,0x04))
	{
		printf("文件不存在或者打开失败! \n");
		return false;
	}

	//设置文件名
	SAFERELEASEARRAY(m_segyFileName)
	int nLen = strlen(segyFileName) + 1;
	m_segyFileName = new char[nLen];
	strcpy(m_segyFileName, segyFileName);
	m_segyFileType = segyFileType ;

	//分析文件
	if(!AnalyseSegyHeader()) 
	{
		//printf("选择的文件不是Seg-Y格式 ! \n");
		return false;
	}

	m_buffer = new char[m_traceBytes];
	m_bSegyFileOpened = true;
	return true;
}

bool CSegyData::CloseSegyFile( )
{
	SAFERELEASEARRAY(m_segyFileName)
		SAFERELEASEARRAY(m_buffer)

		if (m_fileSegy)
		{
			close(m_fileSegy);
			m_fileSegy = -1;
		}

		return true;
}

int       CSegyData::GetTraceBytes()
{
	return m_traceBytes;
}

int		CSegyData::GetSegyTotalTraceNumber()//总道数
{
	return m_dataTraceNo ;
}
int		CSegyData::GetSegySampleNumber()//每道样点数
{
	return m_lineHead.samplenum ;
}
int		CSegyData::GetSegyTimeInterval()//in 毫秒
{
	return m_lineHead.interval/1000;
}
int		CSegyData::GetSegyDataFormat()//	1 = float point(4B), 	2 = fixed point(4B) 3 = fixed point(2B)	4 = fixed point w/gain code(4B)          
{
	return m_lineHead.dataform;
}
int		CSegyData::GetSegyFileType()//0 --- tape file  1----disk file
{
	return m_segyFileType;
}

char*		CSegyData::GetSegyFileName()
{
	return m_segyFileName;
}
ReelHead	CSegyData::GetSegyReelHead() //3200 bytes
{
	return m_reelHead;
}
LineHead	CSegyData::GetSegyLineHead() //400 bytes
{
	return m_lineHead;
}
TraceHead	CSegyData::GetSegyTraceHead() //240 bytes
{
	return m_traceHead;
}

void CSegyData::SwapInt4(char *p)
{
	char c;

	c = p[0]; p[0] = p[3]; p[3] = c;
	c = p[1]; p[1] = p[2]; p[2] = c;
}

void CSegyData::SwapInt4(char *p, int nCount)
{
	char c;
	char *pData = p;
	int i;
	for (i=0; i<nCount; i++)
	{
		c = pData[0]; pData[0] = pData[3]; pData[3] = c;
		c = pData[1]; pData[1] = pData[2]; pData[2] = c;

		pData += 4;
	}
}

unsigned short   CSegyData::SwapInt2(const unsigned short &ll)
{
	char  b;
	unsigned short   aa = ll;
	b = 0xff & (aa>>8) ;
	aa <<= 8;
	aa |= 0x00ff & b ;
	return aa;
}

short CSegyData::SwapInt2(const short &ii)
{
	union{
		short ii;
		unsigned short ui;
	}dd;
	dd.ii = ii;
	char  b;
	b = 0xff & (dd.ui>>8) ;
	dd.ui <<= 8;
	dd.ui |= 0x00ff & b ;
	return dd.ii ;
}

unsigned long CSegyData::SwapInt4(const unsigned long &ll)
{
	unsigned long aa ;
	unsigned short lo,hi;
	lo=unsigned short (0x0000ffff & ll );
	hi =unsigned short (0x0000ffff & (ll>>16) );
	lo = SwapInt2(lo);
	hi = SwapInt2(hi);
	aa = 0 ;
	aa = 0x000ffff & lo ;
	aa <<= 16;
	aa |= 0x0000ffff & hi ;
	return aa;
}
long CSegyData::SwapInt4(const long &ll)
{
	union{
		long ll;
		unsigned long ul;
	}dd;
	dd.ll = ll ;
	unsigned short lo,hi;
	lo=unsigned short(0x0000ffff & dd.ul );
	hi =unsigned short(0x0000ffff &( dd.ul>>16) );
	lo = SwapInt2(lo);
	hi = SwapInt2(hi);
	dd.ul= 0 ;
	dd.ul= 0x000ffff & lo ;
	dd.ul<<= 16;
	dd.ul|= 0x0000ffff & hi ;
	return dd.ll;
}
unsigned int	CSegyData::SwapInt4(const unsigned int &ll)
{
	union{
		int ll;
		unsigned int ul;
	}dd;
	dd.ul = ll ;
	unsigned short lo,hi;
	lo=0xffff & dd.ul ;
	hi = 0xffff & (dd.ul>>16) ;
	lo = SwapInt2(lo);
	hi = SwapInt2(hi);
	dd.ul= 0 ;
	dd.ul= 0x000ffff & lo ;
	dd.ul<<= 16;
	dd.ul|= 0x0000ffff & hi ;
	return dd.ul;
}
int CSegyData::SwapInt4(const int &ll)
{
	union{
		int ll;
		unsigned int ul;
	}dd;
	dd.ll = ll ;
	unsigned short lo,hi;
	lo=0xffff & dd.ul ;
	hi = 0xffff & (dd.ul>>16) ;
	lo = SwapInt2(lo);
	hi = SwapInt2(hi);
	dd.ul= 0 ;
	dd.ul= 0x000ffff & lo ;
	dd.ul<<= 16;
	dd.ul|= 0x0000ffff & hi ;
	return dd.ll;
}

float CSegyData::SwapInt4(const float &ff)
{
	union{
		float ff;
		unsigned long ll;
	}dd;
	dd.ff = ff ;

	unsigned short lo,hi;
	lo=unsigned short (0x0000ffff & dd.ll );
	hi = unsigned short (0x0000ffff & (dd.ll>>16) );
	lo = SwapInt2(lo);
	hi = SwapInt2(hi);
	dd.ll = 0 ;
	dd.ll = 0x000ffff & lo ;
	dd.ll <<= 16;
	dd.ll |= 0x0000ffff & hi ;
	return dd.ff;
}

float CSegyData::ZPOW(const float &x,const int &y)
{
	int i,sign,yy;
	float s;
	yy=y;
	if(yy<0){ yy=-yy;sign=1;}
	else sign=0;
	for(s=1.0,i=0;i<yy;i++) s*=x;
	if(sign) s=float (1.0/s);
	return s;

}

/****************************************************************************
Seg-y to float
****************************************************************************/
float   CSegyData::SegYToFloat( const float &data )
{
	int            index, sign ;
	float          sum = 0. ;
	unsigned long  middle ;
	union {
		float a;
		unsigned long b;
	}c;
	c.a=data;
	middle = c.b >> 31 ;
	if( middle == 1 ) sign = -1 ;
	else  sign = 1 ;
	middle = c.b   >> 24 ;

	index = middle & 127 ;
	index -= 64 ;
	sum =  float ( c.b & 0x00FFFFFF );
	sum /= 16777216.0 ; // pow( 2.0, 24 ) = 16777216
	sum *= float ( pow( 16. , index ) * sign );
	return sum ;
}

unsigned long CSegyData::SegYToFloat( const unsigned long &data )
{
	int            index, sign ;
	float          sum = 0. ;
	unsigned long  middle ;
	union {
		float a;
		unsigned long b;
	}c;
	c.b=data;
	middle = c.b >> 31;
	if( middle == 1 ) sign = -1 ;
	else  sign = 1 ;
	middle = c.b >> 24 ;

	index = middle & 127 ;
	index -= 64 ;
	sum =  float ( c.b & 0x00FFFFFF );
	sum /= 16777216.0 ; // pow( 2.0, 24 ) = 16777216
	sum *= float ( pow( 16. , index ) * sign );
	c.a= sum ;
	return c.b ;

}
long CSegyData::SegYToFloat( const long &data )
{
	int            index, sign ;
	float          sum = 0. ;
	unsigned long  middle ;
	union {
		float f;
		long a;
		unsigned long b;
	}c;
	c.a=data;
	middle = c.b >> 31 ;
	if( middle == 1 ) sign = -1 ;
	else  sign = 1 ;
	middle = c.b   >> 24 ;

	index = middle & 127 ;
	index -= 64 ;
	sum =  float ( c.b & 0x00FFFFFF );
	sum /= 16777216.0 ;      /* pow( 2.0, 24 ) = 16777216 ; */
	sum *= float ( pow( 16. , index ) * sign );
	c.f= sum ;
	return c.a ;
}

unsigned long CSegyData::FloatToSegy(const unsigned long &a)
{
	unsigned long S,M,m ;
	short    int  E, e  ;
	union{
		unsigned long F ;
		float    f ;
	}b;
	b.F=a;

	m=b.F<<9>>9|0x800001;
	e=short(b.F<<1>>24);
	e-=126;

	if((int)e>0){
		E=e+4-e%4;
		M=m>>(4-e%4);
	}
	else{
		E=e-e%4;
		M=m>>(-e%4);
	}
	E/=4;E+=64;

	S&=0x00000000;
	S|=b.F>>31;
	S<<= 7;
	S|=(unsigned long)E;
	S<<=24;
	S|=M;
	return S ;
}

long CSegyData::FloatToSegy(const long &a)
{
	unsigned long S,M,m ;
	short    int  E, e  ;
	union{
		unsigned long F ;
		long	f ;
	}b;
	b.f=a;

	m=b.F<<9>>9|0x800001;
	e=short(b.F<<1>>24);
	e-=126;

	if((int)e>0){
		E=e+4-e%4;
		M=m>>(4-e%4);
	}
	else{
		E=e-e%4;
		M=m>>(-e%4);
	}
	E/=4;E+=64;

	S&=0x00000000;
	S|=b.F>>31;
	S<<= 7;
	S|=(unsigned long)E;
	S<<=24;
	S|=M;
	b.F = S;
	return b.f ;
}

float CSegyData::FloatToSegy(const float &a)
{
	unsigned long S,M,m ;
	short    int  E, e  ;
	union{
		unsigned long F ;
		float    f ;
	}b;
	b.f=a;

	m=b.F<<9>>9|0x800001;
	e=short(b.F<<1>>24);
	e-=126;

	if((int)e>0){
		E=e+4-e%4;
		M=m>>(4-e%4);
	}
	else{
		E=e-e%4;
		M=m>>(-e%4);
	}
	E/=4;E+=64;

	S&=0x00000000;
	S|=b.F>>31;
	S<<= 7;
	S|=(unsigned long)E;
	S<<=24;
	S|=M;
	b.F = S;
	return b.f ;
}

//private member function begin
bool CSegyData::AnalyseSegyHeader( )
{
	if( (m_fileSegy = _open((const char *)m_segyFileName,O_RDONLY|O_BINARY) )==-1  ) 
	{
		printf("Cannot open file : %s !\n",m_segyFileName) ; 
		return false;
	}

	/**************  Reel Head Messages ***************************/
	int i=read(m_fileSegy ,(char *)&m_reelHead,3200);

	if(i<3200) 
	{
		close(m_fileSegy); 
		m_fileSegy = -1;
		printf("File : %s read ReelHead error !\n",(const char *)m_segyFileName) ; 
		return false;
	}

	/*************** Read Line Head Messages ********************/
	i = read( m_fileSegy,(char *)&m_lineHead,400);  

	if(i<400)
	{
		close(m_fileSegy); m_fileSegy = -1;
		printf("File : %s read LineHead error !\n",(const char *)m_segyFileName) ; 
		return false;
	} 

	/*************** Read First Trace Head Messages ********************/
	i = read( m_fileSegy,(char *)&m_traceHead,240); 

	if(i<240)
	{
		close(m_fileSegy);  m_fileSegy = -1;
		printf("File : %s read TraceHead error !\n",(const char *)m_segyFileName) ; 
		return false;
	} 

	//如果是磁带seg-y文件，需要交换各个头结构的字节顺序
	if( m_segyFileType == 0 )  //0 --- tape file  1----disk file
	{
		swapSegyHeadBytes();
	}

	if( (m_lineHead.interval != m_traceHead.SampleInterval) || 
		(m_lineHead.samplenum != m_traceHead.NumberSamples) ||
		((m_lineHead.dataform != 1) && (m_lineHead.dataform !=2) &&
		(m_lineHead.dataform != 3) && (m_lineHead.dataform !=4)  && (m_lineHead.dataform !=5)))
	{
		//AfxMessageBox("Not A Seg-Y File Format ! ");
		//printf("选择的文件不是Seg-Y格式!");
		//return false;
	}

	/*------------------------------------------------*/
	if( m_lineHead.dataform==0||//数据格式缺失时增加, by fei gaolei
		m_lineHead.dataform==1||
		m_lineHead.dataform==2||
		m_lineHead.dataform==4||
		m_lineHead.dataform==5 )
	{
		m_traceBytes = m_lineHead.samplenum*4 + 240;
	}
	else if(m_lineHead.dataform == 3)
	{
		m_traceBytes = m_lineHead.samplenum*2 + 240;
	}

	__int64 filesize;
	filesize = _lseeki64(m_fileSegy,0L,SEEK_END );//SEEK_CUR SEEK_SET
	m_dataTraceNo = ( filesize - 3200 - 400 ) / m_traceBytes;
	return true;
}

void CSegyData::swapSegyHeadBytes()
{
	swapSegyLineHeadBytes( m_lineHead); //400 bytes
	swapSegyTraceHeadBytes(m_traceHead ); //240 bytes
}

void CSegyData::swapSegyLineHeadBytes( LineHead	&lineHead) //400 bytes
{
	/********************************************************************************
	Line Head Messages 
	********************************************************************************/
	lineHead.jobnum = SwapInt4( lineHead.jobnum  ) ;
	/* only one line per reel */
	lineHead.linenum= SwapInt4( lineHead.linenum  ) ;
	lineHead.reelnum= SwapInt4( lineHead.reelnum ) ;
	/* number of data traces per record */
	lineHead.tracenum= SwapInt2( lineHead.tracenum ) ;
	/* number of auxiliary traces per record */
	lineHead.auxtrace= SwapInt2( lineHead.auxtrace ) ;
	/* for this reel of data */
	lineHead.interval= SwapInt2( lineHead.interval ) ; 
	/* for original field recording */
	lineHead.interval0= SwapInt2( lineHead.interval0) ;
	/* for this reel of data */
	lineHead.samplenum= SwapInt2( lineHead.samplenum) ;
	/* for original field recording */
	lineHead.samplenum0= SwapInt2( lineHead.samplenum0 ) ; 
	//	1 = float point(4B), 	2 = fixed point(4B) 3 = fixed point(2B)	4 = fixed point w/gain code(4B)          
	lineHead.dataform= SwapInt2( lineHead.dataform) ;  
	/* expected number of data trace per CDP ensemble */
	lineHead.cdpfold= SwapInt2( lineHead.cdpfold) ;

}// end swapSegyLineHeadBytes


/*****
功能：交换segy格式地震数据道头信息(240 bytes交换字节)
参  数：
traceHead：segyTraceHead类型的结构体
返回值：无
******/
void	CSegyData::swapSegyTraceHeadBytes(TraceHead &traceHead )//240 bytes
{
	traceHead.TraceLsNum	 = SwapInt4( traceHead.TraceLsNum );
	traceHead.TraceRsNum	 = SwapInt4( traceHead.TraceRsNum );
	traceHead.RecordNum	 = SwapInt4( traceHead.RecordNum );
	traceHead.TrNum	 = SwapInt4( traceHead.TrNum );
	traceHead.EnergySourcePointNum	 = SwapInt4( traceHead.EnergySourcePointNum );
	traceHead.CdpNum	 = SwapInt4( traceHead.CdpNum );
	traceHead.TrNumCdp	 = SwapInt4( traceHead.TrNumCdp );
	traceHead.ShotRecOffset	 = SwapInt4( traceHead.ShotRecOffset );
	traceHead.ReceiverElevation	 = SwapInt4( traceHead.ReceiverElevation );
	traceHead.SubfaceElevationAtSource	 = SwapInt4( traceHead.SubfaceElevationAtSource );
	traceHead.SourceDepthBelowSurface	 = SwapInt4( traceHead.SourceDepthBelowSurface );
	traceHead.DatumElevationAtReceiver	 = SwapInt4( traceHead.DatumElevationAtReceiver );
	traceHead.DatumElevationAtSource	 = SwapInt4( traceHead.DatumElevationAtSource );
	traceHead.WaterDepthAtSource	 = SwapInt4( traceHead.WaterDepthAtSource );
	traceHead.WaterDepthAtGroup	 = SwapInt4( traceHead.WaterDepthAtGroup );
	traceHead.SourceX	 = SwapInt4( traceHead.SourceX );
	traceHead.SourceY	 = SwapInt4( traceHead.SourceY );
	traceHead.GroupX	 = SwapInt4( traceHead.GroupX );
	traceHead.GroupY		 = SwapInt4( traceHead.GroupY );	
	traceHead.SourceResidualStaticCor	 = SwapInt4( traceHead.SourceResidualStaticCor );
	traceHead.GroupResidualStaticCor	 = SwapInt4( traceHead.GroupResidualStaticCor );
	traceHead.CmpX	 = SwapInt4( traceHead.CmpX );
	traceHead.CmpY	 = SwapInt4( traceHead.CmpY );
	traceHead.CoordinateProjectType	 = SwapInt4( traceHead.CoordinateProjectType );
	traceHead.WaterDepthAtCdp	 = SwapInt4( traceHead.WaterDepthAtCdp );
	traceHead.DatumTime	 = SwapInt4( traceHead.DatumTime );
	traceHead.DatumVelo  	 = SwapInt4( traceHead.DatumVelo );  


	traceHead.TraceId	 = SwapInt2( traceHead.TraceId );
	traceHead.NumVerStackedTraces	 = SwapInt2( traceHead.NumVerStackedTraces );
	traceHead.NumHorStackedTraces	 = SwapInt2( traceHead.NumHorStackedTraces );
	traceHead.DataUse	 = SwapInt2( traceHead.DataUse );
	traceHead.ElevScale	 = SwapInt2( traceHead.ElevScale );
	traceHead.CoordinateScale	 = SwapInt2( traceHead.CoordinateScale );
	traceHead.CoordinateUnits	 = SwapInt2( traceHead.CoordinateUnits );
	traceHead.WeatheringVelocity	 = SwapInt2( traceHead.WeatheringVelocity );
	traceHead.SubWeatheringVelocity	 = SwapInt2( traceHead.SubWeatheringVelocity );
	traceHead.UpholeTimeAtSource	 = SwapInt2( traceHead.UpholeTimeAtSource );
	traceHead.UpholeTimeAtGroup	 = SwapInt2( traceHead.UpholeTimeAtGroup );
	traceHead.SourceStaticCorrection	 = SwapInt2( traceHead.SourceStaticCorrection );
	traceHead.GroupStaticCorrection	 = SwapInt2( traceHead.GroupStaticCorrection );
	traceHead.TotalStatic	 = SwapInt2( traceHead.TotalStatic );
	traceHead.LagTimeA	 = SwapInt2( traceHead.LagTimeA );
	traceHead.LagTimeB	 = SwapInt2( traceHead.LagTimeB );
	traceHead.RecordingDelay	 = SwapInt2( traceHead.RecordingDelay );
	traceHead.MuteTimeStart	 = SwapInt2( traceHead.MuteTimeStart );
	traceHead.MuteTimeEnd	 = SwapInt2( traceHead.MuteTimeEnd );
	traceHead.NumberSamples	 = SwapInt2( traceHead.NumberSamples );
	traceHead.SampleInterval	 = SwapInt2( traceHead.SampleInterval );
	traceHead.InstrumentGainType	 = SwapInt2( traceHead.InstrumentGainType );
	traceHead.GainConstant	 = SwapInt2( traceHead.GainConstant );
	traceHead.InitialGain	 = SwapInt2( traceHead.InitialGain );
	traceHead.CorrelatedTrace	 = SwapInt2( traceHead.CorrelatedTrace );
	traceHead.SweepFreqAtStart	 = SwapInt2( traceHead.SweepFreqAtStart );
	traceHead.SweepFreqAtEnd	 = SwapInt2( traceHead.SweepFreqAtEnd );
	traceHead.SweepLength	 = SwapInt2( traceHead.SweepLength );
	traceHead.SweepType	 = SwapInt2( traceHead.SweepType );
	traceHead.SweepStartTaperLength	 = SwapInt2( traceHead.SweepStartTaperLength );
	traceHead.SweepEndTaperLength	 = SwapInt2( traceHead.SweepEndTaperLength );
	traceHead.TaperType	 = SwapInt2( traceHead.TaperType );
	traceHead.AliasFilterFreq	 = SwapInt2( traceHead.AliasFilterFreq );
	traceHead.AliasFilterSlope	 = SwapInt2( traceHead.AliasFilterSlope );
	traceHead.NotchFilterFreq	 = SwapInt2( traceHead.NotchFilterFreq );
	traceHead.NotchFilterSlope	 = SwapInt2( traceHead.NotchFilterSlope );
	traceHead.LowCutFreq	 = SwapInt2( traceHead.LowCutFreq );
	traceHead.HighCutFreq	 = SwapInt2( traceHead.HighCutFreq );
	traceHead.LowCutSlope	 = SwapInt2( traceHead.LowCutSlope );
	traceHead.HighCutSlope	 = SwapInt2( traceHead.HighCutSlope );
	traceHead.Year	 = SwapInt2( traceHead.Year );
	traceHead.Day	 = SwapInt2( traceHead.Day );
	traceHead.Hour	 = SwapInt2( traceHead.Hour );
	traceHead.Minute	 = SwapInt2( traceHead.Minute );
	traceHead.Second	 = SwapInt2( traceHead.Second );
	traceHead.TimeCode	 = SwapInt2( traceHead.TimeCode );
	traceHead.TrWeightFactor	 = SwapInt2( traceHead.TrWeightFactor );
	traceHead.GroupNumberOfRollSwitch	 = SwapInt2( traceHead.GroupNumberOfRollSwitch );
	traceHead.GroupNumOf1stTrcInOrig	 = SwapInt2( traceHead.GroupNumOf1stTrcInOrig );
	traceHead.GroupNumOfLastTrcInOrig	 = SwapInt2( traceHead.GroupNumOfLastTrcInOrig );
	traceHead.GapSize	 = SwapInt2( traceHead.GapSize );
	traceHead.TaperOvertravely	 = SwapInt2( traceHead.TaperOvertravely );
	traceHead.DataTrnum	 = SwapInt2( traceHead.DataTrnum );//若DataTrnum为short型则修改回SwapInt2
	traceHead.CmpDatumStaticCor	 = SwapInt2( traceHead.CmpDatumStaticCor );//修改的地方
	traceHead.GroupStaNum	 = SwapInt2( traceHead.GroupStaNum );
	traceHead.Linenum3dBack	 = SwapInt2( traceHead.Linenum3dBack );
	traceHead.FieldFileNum	 = SwapInt2( traceHead.FieldFileNum );
	traceHead.SourceStaNum	 = SwapInt2( traceHead.SourceStaNum );
	traceHead.CmpElevation	 = SwapInt2( traceHead.CmpElevation );
	traceHead.CmpStaNum	 = SwapInt2( traceHead.CmpStaNum );
	traceHead.TMin	 = SwapInt2( traceHead.TMin );
	traceHead.TMax	 = SwapInt2( traceHead.TMax );
	traceHead.Dimension	 = SwapInt2( traceHead.Dimension );
}




