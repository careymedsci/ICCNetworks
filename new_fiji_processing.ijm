# ============================================================================
# Copyright (c) 2025 Thomas Broggini and Liu Xiao. All rights reserved.
#
# This code was jointly written by:
# Thomas Broggini (Frankfurt University, Germany)
# Liu Xiao (Xiangyang First Hospital, China)
#
# The way how the mathematical models and methodologies implemented in this code have been
# individually customized and are not intended for generic use. 
#
# For permission requests or inquiries, please contact:
# Thomas Broggini  : broggini@med.uni-frankfurt.de
# Liu Xiao         : careyneurosurgery@gmail.com  /  carey-lau@foxmail.com
#
# Unauthorized use will be considered a violation of intellectual property rights.
# ============================================================================




   //作者：刘晓；Thomas Broggini; 法兰克福大学医院；湖北医药学院襄阳市第一人民医院
   //By Xiao Liu and Thomas Broggini； University Clinic Frankfurt； Xiangyang No.1 people's Hospital (Xiao Liu)
   //2024.10.09
   //Also an matlab version coded by Xiao Liu; Xiao.Liu@stud.uni-frankfurt.de


   //打开要处理的视频
   //opens folder where the denoised video sits in

#@ File (label = "Input directory", style = "directory") inputdir
   list = getFileList(inputdir);
   list = Array.sort(list);
   output = inputdir + "_roi_centre";
   outputrois = inputdir + "_rois_shape";
   File.makeDirectory(output);
   File.makeDirectory(outputrois);
   
   //列举文件夹内所有待处理视频
   //Lists all videos
   for (i = 0; i < list.length; i++) { 
	xx = i +1;
	print("Time Series " + xx + " of " + list.length);

    //清空ROI管理器
	//clears Roi manager
	nROIs = roiManager("count");
	if (nROIs > 0)
	{
		roiManager("Deselect");
		roiManager("Delete");
	}

    //将所有帧中最大像素映射到一张图片上
	//opens file and creates Max_z_projection
	open(list[i]);
	filename = list[i];
	run("8-bit");
	
	//run("Brightness/Contrast...");
    run("Enhance Contrast", "saturated=0.35");
    run("Apply LUT", "stack");
    run("Median 3D...", "x=2 y=2 z=2");

    run("Z Project...", "projection=[Max Intensity]");

    //run("Brightness/Contrast...");
    run("Enhance Contrast", "saturated=0.35");
    run("Apply LUT");
    
     //转为灰度图像
	//Converts to grayscale
	run("Conversions...", "scale");
	
	//保存最大像素图片
    //saves max_z_project
	run("Duplicate...", " ");
	index = lastIndexOf(filename, "."); 
  	if (index!=-1) filename2 = substring(filename, 0, index);
	saveAs("tif", output + "/" + filename2 + "_max_z_proj.tif");	
	close();
	
	waitForUser("please check the picture，press OK to contiune");
    print("press OK to contiune...");

    //二掩膜
    //Create biniary Photo
	run("Auto Local Threshold", "method=Otsu radius=60 parameter_1=0 parameter_2=0 white");
	setOption("BlackBackground", true);
	
	//细分ROI区
	//Watershed
	run("Watershed", "stack");
	
	//添加识别的ROI到管理器
	//adds all rois to manager
	run("Analyze Particles...", "size=50-800 add");
	close();
	
	//展示所有ROI
	// Shows all ROIs
	roiManager("Show All");
	
	waitForUser("check the image, press ok to contiune");
    // 执行后续操作
    print("countiue...");
	
	//获取ROI数量
	//Lists all Rois
    roiCount = roiManager("count");
    count = 0;
    //循环保存坐标，遍历所有现有ROI
    //caculates all Rois
    for (i = 0; i < roiCount; i++) {
    roiManager("Select", i);

    //生成文件名
    //Names the file of result
    filename3 = outputrois + "/" + "roi_coordinates_" + count + ".txt"; 
    //保存ROI坐标
    //save all detailed rois with orignal shape
    run("Save XY Coordinates...", "save=[" + filename3 + "]");
    count++;
    }
    
	//展示所有ROI
	// Shows all ROIs
	//close();
	roiManager("Show All");
	
    //重绘所有感兴趣区，以便进一步精确
	//scales all ROIS with factor 0.8
	for (j = 0; j < roiManager("count"); j++) {
		roiManager("select", j);
		run("Scale... ", "x=0.8 y=0.8 centered");
		roiManager("update")
		}
    
    //记录所有单点坐标
	//Measures and saves location of each ROIs（central）
	count=roiManager("count"); 
	array=newArray(count); 
	for(k=0; k<count;k++) { 
        array[k] = k; 
	} 
	roiManager("Select", array);
	run("Select All");
	run("Set Measurements...", "center redirect=None decimal=0");
	roiManager("Multi Measure");
  	selectWindow("Results");
	saveAs("Results", output + "/" + filename2 + "_coordintes.csv");
	close("Results");

    //测量所有ROI的平均值
	//Measures mean intensities in each ROIs over time
	roiManager("Show None");
	roiManager("Show All");
	run("Set Measurements...", "mean redirect=None decimal=0");
	roiManager("Multi Measure");
	
    //保存ROI的平均像素值
	//saves peixl values results as .txt file in output folder and closes results
  	selectWindow("Results");
	saveAs("Results", output + "/" + filename2 + "_mean_values.csv");
	close("Results");

	close("*");


}