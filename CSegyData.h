// SegyData.h: interface for the CSegyData class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CSEGYDATA_H
#define CSEGYDATA_H
#if _MSC_VER > 1000
#pragma once
#endif 

#include <stdlib.h>
#include <stdio.h>
#include <io.h>
#include <string.h>
#include <cstring>
#ifndef SAFERELEASE
#define SAFERELEASE(X) {if(X){delete X; X=NULL;}}
#define SAFERELEASEARRAY(X) {if(X){delete[] X; X=NULL;}}
#endif

#define  Float32	1	/* Beyong to Binary */
#define  Int32		2	/* Beyong to Binary */
#define  Int16		3	/* Beyong to Binary */
#define  Int32W		4	/* Beyong to Binary */

/****************** tape head describe *************************************/
struct ReelHead
{
	struct{
		char card[5];          /*= "C01 "*/
		char client0[8];       /*= "CLIENT "*/
		char client[20];
		char company0[10];     /*= "COMPANY "*/
		char company[20];
		char crewno0[10];      /*= "CROW NO "*/
		char crewno[7];
	}    card1;

	struct{
		char card[5];      /*="C02 "*/
		char line0[6];      /*="LINE "*/
		char line[20];
		char area0[6];      /*="AREA "*/
		char area[20];
		char mapid0[8];      /*="MAP ID "*/
		char mapid[15];
	}    card2;

	struct{
		char card[5];      /*="C03 "*/
		char reelno0[9];      /*="REEL NO "*/
		char reelno[20];
		char cc[25];
		char observer0[11];      /*="OBSERVER "*/
		char observer[10];
	}    card3;

	struct{
		char card[5];      /*="C04 "*/
		char cc[75];
	}    card4;

	struct{
		char card[5];      /*="C05 "*/
		char c1[20];      /*="DATA TRACES/RECORD "*/
		char traceperrecord[5];
		char auxtrace0[11];      /*="AUX TRACE "*/
		char auxtrace[9];
		char sumfold0[10];      /*="SUM FOLD "*/
		char sumfold[5];
		char cdpfold0[10];      /*="CDP FOLD "*/
		char cdpfold[5];
	}    card5;
	struct{
		char card[5];      /*="C06 "*/
		char c1[18];      /*="SAMPLES INTERVAL "*/
		char interval[10];
		char c2[14];      /*="SAMPES/TRACE "*/
		char samplesnum[10];
		char c3[14];      /*="BYTES/SAMPLE "*/
		char size[9];
	}    card6;
	struct{
		char card[5];      /*="C07 "*/
		char c1[18];      /*="RECORDING FORMAT "*/
		char format0[7];
		char c2[18];      /*="FORMAT THIS REEL "*/
		char format[7];
		char c3[20];      /*="MEASUREMENT SYSTEM "*/
		char measure[5];
	}    card7;
	struct{
		char card[5];      /*="C08 "*/
		char c1[23];      /*="SAMPLE CODE  FLOATING "*/
		char samplecode[10];
		char c2[42];
	}    card8;
	struct{
		char card[5];      /*="C09 "*/
		char c1[11];      /*="GAIN TYPE "*/
		char gaintype[20];
		char c2[12];      /*="GAIN LEVEL "*/
		char gainlevel[32];
	}card9;
	struct{
		char card[5];      /*="C10 "*/
		char c[75];
	}    card10;
	struct{
		char card[5];      /*="C11 "*/
		char c[75];
	}    card11;
	struct{
		char card[5];      /*="C12 "*/
		char c[75];
	}    card12;
	struct{
		char card[5];      /*="C13 "*/
		char c[75];
	}    card13;
	struct{
		char card[5];      /*="C14 "*/
		char c[75];
	}    card14;
	struct{
		char card[5];      /*="C15 "*/
		char c[75];
	}    card15;
	struct{
		char card[5];      /*="C16 "*/
		char c[75];
	}    card16;
	struct{
		char card[5];      /*="C17 "*/
		char c[75];
	}    card17;
	struct{
		char card[5];      /*="C18 "*/
		char c1[20];      /*="TRACES SORTED BY : "*/
		char sort[55];
	}    card18;
	struct{
		char card[5];      /*="C19 "*/
		char c[75];      /*="AMPLITUDE RECOVERY "*/
	}    card19;
	struct{
		char card[5];      /*="C20 "*/
		char c1[16];      /*="MAP PROJECTION "*/
		char project[14];
		char c2[10];      /*="ZONE ID "*/
		char zone[10];
		char c3[18];      /*="COORDINATE UNITS "*/
		char units[7];
	}card20;
	struct{
		char card[5];      /*="C21 "*/
		char c[75];
	}card21;
	struct{
		char card[5];      /*="C22 "*/
		char c[75];
	}card22;
	struct{
		char card[5];      /*="C23 "*/
		char c[75];
	}card23;
	struct{
		char card[5];      /*="C24 "*/
		char c[75];
	}card24;
	struct{
		char card[5];      /*="C25 "*/
		char c[75];
	}card25;
	struct{
		char card[5];      /*="C26 "*/
		char c[75];
	}card26;
	struct{
		char card[5];      /*="C27 "*/
		char c[75];
	}card27;
	struct{
		char card[5];      /*="C28 "*/
		char c[75];
	}card28;
	struct{
		char card[5];      /*="C29 "*/
		char c[75];
	}card29;
	struct{
		char card[5];      /*="C30 "*/
		char c[75];
	}card30;
	struct{
		char card[5];      /*="C31 "*/
		char c[75];
	}card31;
	struct{
		char card[5];      /*="C32 "*/
		char c[75];
	}card32;
	struct{
		char card[5];      /*="C33 "*/
		char c[75];
	}card33;
	struct{
		char card[5];      /*="C34 "*/
		char c[75];
	}card34;
	struct{
		char card[5];      /*="C35 "*/
		char c[75];
	}card35;
	struct{
		char card[5];      /*="C36 "*/
		char c[75];
	}card36;
	struct{
		char card[5];      /*="C37 "*/
		char c[75];
	}card37;
	struct{
		char card[5];      /*="C38 "*/
		char c[75];
	}card38;
	struct{
		char card[5];      /*="C39 "*/
		char c[75];
	}card39;
	struct{
		char card[5];      /*="C40 "*/
		char c[75];      /*="END EBCDIC "*/
	}card40;
};                  



/******************************Line Head describe ***************************/
struct  LineHead
{
	long int jobnum;
	long int linenum;   /* only one line per reel */
	long int reelnum;
	short      tracenum;  /* number of data traces per record */
	short      auxtrace;  /* number of auxiliary traces per record */
	short      interval; /* for this reel of data */
	short      interval0;  /* for original field recording */
	short      samplenum;/* for this reel of data */
	short      samplenum0; /* for original field recording */
	short      dataform;  /* 1 = float point(4B), 2 = fixed point(4B) */
	/* 3 = fixed point(2B)                      */
	/* 4 = fixed point w/gain code(4B)          */
	short      cdpfold;   /* expected number of data trace per CDP ensemble */
	char     rest[372];
};

/***************************** Trace Head describe ***************************/
struct TraceHead
{
	int   TraceLsNum;                      /* 1:     1-4    一条测线中的道顺序号。如果一条测线有若干卷磁带，顺序号连续递增*/
	int   TraceRsNum;                     /* 2:     5-8    在本卷磁带中的道顺序号。每卷带的道顺序*/
	int   RecordNum;                      /* 3:     9-12   号从1 开始。原始的野外记录号 */
	int   TrNum;                          /* 4:    13-16   在原始野外记录中的道号*/
	int   EnergySourcePointNum;           /* 5:    17-20   震源点号（在同一个地面点有多于一个记录时使用）*/
	int   CdpNum;                         /* 6:    21-24   CMP 号*/
	int   TrNumCdp;                       /* 7:    25-28   在CMP 道集中的道号（在每个CMP 道集中道号从1 开始） */
	short  TraceId;                           /* 8:    29-30   道识别码：
										   1=地震数据 4=时断 7=记录
										   2=死道 5= 井口时间 8=水断
										   3=DUMMY 6=扫描道 9~N=选择使用（N=32767）
										   */
	short  NumVerStackedTraces;            /* 9:    31-32  产生这一道的垂直叠加道数(1 是一道；2 是两道相加；…) */
	short  NumHorStackedTraces;            /* 10:   33-34  产生这一道的水平叠加道数(1 是一道；2 是两道相加；…)*/
	short  DataUse;                        /* 11:   35-36  数据类型: 1=生产 2=试验 */
	int   ShotRecOffset;                  /* 12:   37-40  从炮点到接收点的距离（如果是相反向激发为负值） */
	int   ReceiverElevation;              /* 13:   41-44  接收点高程。高于海平面的高程为正，低于海平面为负*/
	int   SubfaceElevationAtSource;       /* 14:   45-48  炮点的地面高程*/
	int   SourceDepthBelowSurface;        /* 15:   49-52  炮点低于地面的深度（正数）*/
	int   DatumElevationAtReceiver;       /* 16:   53-56  接收点的基准面高程*/
	int   DatumElevationAtSource;         /* 17:   57-60  炮点的基准面高程*/
	int   WaterDepthAtSource;             /* 18:   61-64  炮点的水深*/
	int   WaterDepthAtGroup;              /* 19:   65-68  接收点的水深*/
	short  ElevScale;                      /* 20:   69-70  对41～68 字节中的所有高程和深度应用了此因子给出真值。
										   比例因子＝1，±10，±100，±l000 或者±l0000。如果为正，乘以因子；
										   如果为负，则除以因子*/
	short  CoordinateScale;                /* 21:   71-72  如果为负，则除以因子。
										   对73～88 字节中的所有坐标应用了此因子给出真值。
										   比例因子＝1，±10，±100，±l000 或者±l0000。
										   如果为正，乘以因子；如果为负，则除以因子*/
	int   SourceX;                        /* 22:   73-76  炮点坐标----X 
										  如果坐标单位是弧度或秒，X 值代表经度，Y 值代表纬度。
										  正值代表格林威治子午线东或者赤道北的秒数。
										  负值则为西或者南的秒数*/
	int   SourceY;                        /* 23:   76-80  炮点坐标----Y*/
	int   GroupX;                         /* 24:   81-84  接收点坐标----X*/
	int   GroupY;                         /* 25:   85-88  接收点坐标----Y*/
	short  CoordinateUnits;                /* 26:   89-90  坐标单位：1＝长度（米或者英尺）；2＝弧度或秒*/
	short  WeatheringVelocity;             /* 27:   91-92  风化层速度*/
	short  SubWeatheringVelocity;          /* 28:   93-94  风化层下的速度*/
	short  UpholeTimeAtSource;             /* 29:   95-96  震源处的井口时间*/
	short  UpholeTimeAtGroup;              /* 30:   97-98  接收点处的井口时间*/
	short  SourceStaticCorrection;         /* 31:  99-100  炮点的静校正*/
	short  GroupStaticCorrection;          /* 32: 101-102  接收点的静校正*/
	short  TotalStatic;                    /* 33: 103-104  应用的总静校正量(如果没有应用静校正为零)*/
	short  LagTimeA;                       /* 34: 105-106  延迟时间-A，以ms 表示。240 字节的道标识的结束和时间信号之间的时间。
										   如果时间信号出现在道头结束之前为正。
										   如果时间信号出现在道头结束之后为负。
										   时间信号就是起始脉冲，它记录在辅助道上或者由记录系统指定*/
	short  LagTimeB;                      /* 35: 107-108  时间延迟-B，以ms 表示。为时间信号和能量起爆之间的时间。可正可负*/
	short  RecordingDelay;                /* 36: 109-110  时间延迟时间，以ms 表示。
										  能量源的起爆时间和开始记录数据样点之间的时间
										  （深水时，数据记录不从时间零开始。）*/
	short  MuteTimeStart;                 /* 37: 111-112  起始切除时间*/
	short  MuteTimeEnd;                   /* 38: 113-114  结束切除时间*/
	short  NumberSamples;                 /* 39: 115-116  本道的采样点数*/
	short  SampleInterval;                /* 40: 117-118  本道的采样间隔、以 us 表示*/
	short  InstrumentGainType;            /* 41: 119-120  野外仪器的增益类型：
										  1=固定增益 2=二进制增益
										  3=浮点增益 4~N=选择使用
										  */
	short  GainConstant;                  /* 42: 121-122  仪器增益常数*/
	short  InitialGain;                   /* 43: 123-124  仪器起始增益(dB)*/
	short  CorrelatedTrace;               /* 44: 125-126  相关码 1=没有相关 2=相关*/
	short  SweepFreqAtStart;              /* 45: 127-128  起始扫描频率*/
	short  SweepFreqAtEnd;                /* 46: 129-130  结束扫描频率*/
	short  SweepLength;                   /* 47: 131-132  扫描长度，以ms 表示*/
	short  SweepType;                     /* 48: 133-134  扫描类型： 1=线性 2=抛物线 3=指数 4=其他 */
	short  SweepStartTaperLength;         /* 49: 135-136  扫描道起始斜坡长度，以ms 表示 */
	short  SweepEndTaperLength;           /* 50: 137-138  扫描道终了斜坡长度，以ms 表示  */
	short  TaperType;                     /* 51: 139-140  斜坡类型： 1=线性 2=COSP2       3=其他  */
	short  AliasFilterFreq;               /* 52: 141-142  滤假频的频率（如果使用）*/
	short  AliasFilterSlope;              /* 53: 143-144  滤假频的陡度*/
	short  NotchFilterFreq;               /* 54: 145-146  陷波陡率（如果使用）*/
	short  NotchFilterSlope;              /* 55: 147-148  陷波陡度*/
	short  LowCutFreq;                    /* 56: 149-150  低截频率（如果使用）*/
	short  HighCutFreq;                   /* 57: 151-152  高截频率（如果使用）*/
	short  LowCutSlope;                   /* 58: 153-154  低截频率陡度*/
	short  HighCutSlope;                  /* 59: 155-156  高截频率陡度*/
	short  Year;                          /* 60: 157-158 =0  数据记录的年*/
	short  Day;                           /* 61: 159-160 =0  日*/
	short  Hour;                          /* 62: 161-162 =0  小时（24 时制）*/
	short  Minute;                        /* 63: 163-164 =0  分*/
	short  Second;                        /* 64: 165-166 =0  秒*/
	short  TimeCode;                      /* 65: 167-168     时间代码： 1=当地时间 2=格林威治时间 3=其他*/
	short  TrWeightFactor;                /* 66: 169-170     道加权因子。（最小有效位定义为2**（―N），N=0，1，2，…，32767）*/
	short  GroupNumberOfRollSwitch;       /* 67: 171-172 =0  覆盖开关位置1 的检波器道号*/
	short  GroupNumOf1stTrcInOrig;        /* 68: 173-174 =0  在原始野外记录中道号1 的检波器号*/
	short  GroupNumOfLastTrcInOrig;       /* 69: 175-176 =0  在原始野外记录中道号最后一道的检波器号*/
	short  GapSize;                       /* 70: 177-178     缺口大小（滚动的总道数）   */
	short  TaperOvertravely;              /* 71: 179-180     在测线的开始或者结束处的斜坡位置：*/
	 //int  DataTrnum;									  //1=在后面 2=在前面																	*/
	short  DataTrnum;                     /* 72: 181-182    total traces of data */
	char   Tmp1[2];                       /* 73: 183-184    not used */
	int   SourceResidualStaticCor;       /* 74: 185-188    */
	int   GroupResidualStaticCor;        /* 75: 189-192    */
	short  CmpDatumStaticCor;             /* 76: 193-194    */
	short  GroupStaNum;                   /* 77: 195-196    */               
	short  Linenum3dBack;                 /* 78: 197-198    */
	char   Tmp2[2];					      /* 79: 199-200    not used */
	int   CmpX;                          /* 80: 201-204    */
	int   CmpY;                          /* 81: 205-208    */
	int   CoordinateProjectType;         /* 82: 209-212    */
	int   WaterDepthAtCdp;               /* 83: 213-216    */
	short  FieldFileNum;                  /* 84: 217-218    */
	short  SourceStaNum;                  /* 85: 219-220    */
	short  CmpElevation;                  /* 86: 221-222    */
	short  CmpStaNum;                     /* 87: 223-224    */
	short  TMin;                          /* 88: 225-226    */
	short  TMax;                          /* 89: 227-228    */
	int   DatumTime;                     /* 90: 229-232    */
	int   DatumVelo;                     /* 91: 233-236    */
	char   Tmp3[2];                       /* 92: 237-238    not used */
	short  Dimension;                     /* 93: 239-240:   0 for 2D and 1 for 3D*/
};

class CSegyData  
{
public:
	CSegyData();
	CSegyData(char* segyFileName,int segyFileType);
	virtual ~CSegyData();

public:

	char*		m_segyFileName;
    
	ReelHead	m_reelHead;//3200 bytes
	LineHead	m_lineHead;//400 bytes
	TraceHead	m_traceHead;//240 bytes

	int			m_dataTraceNo; //seg-y file total data trace number

	__int64			m_traceBytes; //bytes per trace
	int			m_segyFileType;//0 --- tape file  1----disk file
	void		swapSegyTraceHeadBytes(TraceHead &trackHead );//240 bytes
private:
	char *			m_buffer;
	bool            m_bSegyFileOpened;

	int  m_fileSegy;

public:
	int     GetSegyTotalTraceNumber();//总道数
	int     GetSegySampleNumber();//每道样点数
	int     GetSegyTimeInterval();//in 毫秒
	int     GetSegyDataFormat();//	1 = float point(4B), 	2 = fixed point(4B) 3 = fixed point(2B)	4 = fixed point w/gain code(4B)          
	int     GetSegyFileType();//0 --- tape file  1----disk file
	char*   GetSegyFileName();

	int      GetTraceBytes();

	ReelHead     GetSegyReelHead();//3200 bytes
	LineHead     GetSegyLineHead();//400 bytes
	TraceHead    GetSegyTraceHead();//240 bytes

	bool     OpenSegyFile(const char* segyFileName,const int segyFileType);
	bool     CloseSegyFile( void );

	//根据输入道序号，读取一个地震道
	//参数 int nTraceIndex       : 输入道序号，从0开始计数( 范围： 0 <--> GetSegyTotalTraceNumber() - 1 )
	//参数 TraceHead *pTraceHead : 返回道头指针
	//参数 float *pTraceData     : 返回地震数据数组指针（长度为GetSegySampleNumber()）
	//返回值：是否读取成功，通常，如果输入道序号范围正确，总是会返回true
	bool     ReadOneTrace(__int64 nTraceIndex, TraceHead *pTraceHead, float *pTraceData);

	//根据输入道序号，读取一个地震道
	//参数 int nTraceIndex         : 输入道序号，从0开始计数( 范围： 0 <--> GetSegyTotalTraceNumber() - 1 )
	//参数 float *pTraceData     : 返回地震数据数组指针（长度为GetSegySampleNumber()）
	//返回值：是否读取成功，通常，如果输入道序号范围正确，总是会返回true
	bool     ReadOneTrace(__int64 nTraceIndex, float *pTraceData);

	//根据输入道序号，读取一个地震道
	//参数 int nTraceIndex         : 输入道序号，从0开始计数( 范围： 0 <--> GetSegyTotalTraceNumber() - 1 )
	//参数 TraceHead *pTraceHead : 返回道头指针
	//返回值：是否读取成功，通常，如果输入道序号范围正确，总是会返回true
	bool     ReadOneTraceHeader(__int64 nTraceIndex, TraceHead *pTraceHead);
   /* bool     ReadTraceHeadPartInfo()*/
	bool	SplitSegyFile(const char *newSegyFileName ,
		const __int64 		m_jumpSamp,
		const __int64 		m_jumpTrace,
		const    short 		m_showSamp,
		const    int		m_showTrace	);
	bool	SegyToBinary(char *  w_buffer,//it's length must long enough ;
		const int &		m_jumpSamp,
		const int &		m_jumpTrace,
		const int &		m_showSamp,
		const int &		m_showTrace	);
	bool SegyToBinary(
		const __int64 		m_jumpSamp,
		const __int64		m_jumpTrace,
		const int &		m_showSamp,
		const int &		m_showTrace,
		const char * newBinFileName);
	bool	BinaryToSegy(const char* w_BINFileName,
		const int &		w_binDataType,  //1  --- 32bit float, 2----32bit integer 3----16bit integer
		const int &		w_segyFileType, //0 --- tape file  1----disk file
		const char *newSegyFileName,
		const int &		w_sampNum,
		const int &		w_startCDP,
		const int &		w_traceNum,
		const int &		w_timer);
	//bool WriteSegyFileHead(ReelHead & referReelHead, LineHead &referLineHead,);
	//bool Writer1TraceData2SegyFileEnd(int fdo, char * pCurrentTraceBuffer, long & nTraceIndex, int outTraceBytes, int m_segyFileType);
 //   
public:

	unsigned long	SegYToFloat( const unsigned long &data );
	long			SegYToFloat( const long &data );
	float			SegYToFloat( const float &data );

	unsigned long	FloatToSegy(const unsigned long &a);
	long			FloatToSegy(const long &a);
	float			FloatToSegy(const float &a);

	float			    SwapInt4(const float &ff);
	int				    SwapInt4(const int &ll);
	unsigned int	SwapInt4(const unsigned int &ll);
	long			SwapInt4(const long &ll);
	unsigned long	SwapInt4(const unsigned long &ll);
	short			SwapInt2(const short &ii);
	unsigned short	SwapInt2(const unsigned short &ll);
	void            SwapInt4(char *p);
	void            SwapInt4(char *p, int nCount);

public:

	float			ZPOW(const float &x,const int &y);
	bool			AnalyseSegyHeader( );
	void			swapSegyHeadBytes( );
	void	        swapSegyLineHeadBytes( LineHead	&lineHead);//400 bytes
};


#endif 
