// Cryo-Confocal Collaborator Export Macro
// Fiji / ImageJ 1.54p (requires Bio-Formats)
//
// Downstream export after target screening. Expects in one folder:
//   sample.czi + sample_targets.zip  (from the screening macro)
//
// Workflow (sequential):
//   1. Max-projection PNGs — confocal and reflection, target marks burned in
//   2. Interactive 2000×2000 px region crops (confocal + reflection, same px coords, marks visible)
//   3. Interactive 350×350 px zooms on the target Z-slice (confocal, no marks, 50 µm scale bar), per target

macro "Cryo Confocal Collaborator Export" {
	requires("1.52n");

	inputDir = getDirectory("Select folder containing CZI and targets.zip");
	if (inputDir == "") exit("No input folder selected.");
	if (!endsWith(inputDir, File.separator)) inputDir = inputDir + File.separator;

	outputDir = getDirectory("Select output folder");
	if (outputDir == "") exit("No output folder selected.");
	if (!endsWith(outputDir, File.separator)) outputDir = outputDir + File.separator;

	cziPath = findCziFile(inputDir);
	if (cziPath == "") exit("No .czi file found in:\n" + inputDir);

	baseName = getBaseName(cziPath);
	zipPath = findTargetsZip(inputDir, baseName);
	if (zipPath == "") exit("No targets zip found for:\n" + baseName + "\nin:\n" + inputDir);

	run("Bio-Formats Importer", "open=[" + cziPath + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	stackID = getImageID();
	getDimensions(impW, impH, numChannels, numZ, numFrames);

	if (numChannels != 2 && numChannels != 3) {
		exit("Expected 2 channels (blue + reflection) or 3 (blue + green + reflection).\nFound: " + numChannels);
	}

	hasGreen = (numChannels == 3);
	reflCh = 2;
	if (hasGreen) reflCh = 3;

	// 2-channel: C1 = blue (Hoechst), C2 = reflection. No composite view.
	selectImage(stackID);
	if (hasGreen) {
		run("Make Composite");
	} else {
		if (is("composite")) run("Make Composite");
	}

	roiManager("reset");
	roiManager("Associate", "true");
	roiManager("Open", zipPath);
	selectImage(stackID);

	// ---- Load ROI data (must stay in main block; function-local assignment breaks globals) ----
	roiCount = roiManager("count");
	if (roiCount == 0) exit("No ROIs found in the targets zip.");

	roiIsLabel = newArray(roiCount);
	roiIsTarget = newArray(roiCount);
	roiX = newArray(roiCount);
	roiY = newArray(roiCount);
	roiW = newArray(roiCount);
	roiH = newArray(roiCount);
	roiText = newArray(roiCount);
	targetCount = 0;
	targetCenterX = newArray(500);
	targetCenterY = newArray(500);
	targetZ = newArray(500);

	for (i = 0; i < roiCount; i++) {
		roiManager("select", i);
		roiName = Roi.getName();
		getSelectionBounds(rx, ry, rw, rh);
		roiIsLabel[i] = 0;
		if (indexOf(roiName, "Label_") == 0) roiIsLabel[i] = 1;
		roiIsTarget[i] = 0;
		if (roiIsLabel[i] == 0 && indexOf(roiName, "Target_") == 0) roiIsTarget[i] = 1;
		roiX[i] = rx;
		roiY[i] = ry;
		roiW[i] = rw;
		roiH[i] = rh;
		roiText[i] = "";
		if (roiIsLabel[i] == 1) roiText[i] = getTextRoiText();

		if (roiIsTarget[i] == 1) {
			targetCenterX[targetCount] = rx + rw / 2;
			targetCenterY[targetCount] = ry + rh / 2;
			targetZ[targetCount] = parseTargetDrawnZ(roiName);
			targetCount++;
		}
	}
	run("Select None");

	// Geometry is stored in arrays; remove screening ROIs from stack and manager
	clearAllRoisFromImage(stackID);

	if (targetCount == 0) {
		exit("No Target_* oval ROIs found in:\n" + zipPath + "\n\nRe-run the screening macro or check the zip.");
	}

	// Fixed square sizes in pixels (8045×8045 images; no calibration needed)
	side350Px = 2000;
	side80Px = 350;

	if (hasGreen) channelDesc = " (blue, green, reflection)";
	else channelDesc = " (blue, reflection)";

	showMessage(
		"Collaborator export started.\n\n" +
		"Sample: " + baseName + "\n" +
		"Channels: " + numChannels + channelDesc + "\n" +
		"Targets loaded: " + targetCount + "\n" +
		"Square ROI sizes: " + side350Px + " px (region crop), " + side80Px + " px (target zoom)\n\n" +
		"Do not click Cancel on prompts (that stops the macro)."
	);

	// ---- 1. Confocal max projection + marks ----
	confocalMaxID = buildConfocalMaxProjection(stackID, hasGreen, numZ);
	selectImage(confocalMaxID);
	confocalMaxID = burnStoredRoisOnImage("CONFOCAL_MAX", roiCount);
	saveAs("PNG", outputDir + baseName + "_max_confocal.png");

	// ---- 2. Reflection max projection + B&C + marks (grayscale, single channel) ----
	reflectionMaxID = buildReflectionMaxProjection(stackID, reflCh, numZ);
	run("Brightness/Contrast...");
	waitForUser("Reflection max projection",
		"Adjust Brightness/Contrast for the reflection channel.\n" +
		"Close the B&C window when done, then click OK.\n" +
		"These settings will be used for all reflection exports.");
	selectImage(reflectionMaxID);
	reflectionMaxID = burnStoredRoisOnImage("REFLECTION_MAX", roiCount);
	saveAs("PNG", outputDir + baseName + "_max_reflection.png");

	// ---- 3. Interactive 2000×2000 px region crops ----
	// Position the square on the confocal overview (with marks), read its px coords,
	// then Image>Duplicate the same rectangle on confocal and reflection.
	run("ROI Manager...");
	roiManager("Associate", "false");
	cropCsvPath = outputDir + baseName + "_intermediate_crops.csv";
	cropCsvReady = 0;

	crop350Num = 1;
	while (true) {
		roiManager("reset");
		selectImage(confocalMaxID);
		setTool("rectangle");
		positionSquareRoi(side350Px);
		roiManager("Add");
		setOption("Show All", true);
		roiManager("Show All");
		roiManager("Select", 0);

		waitForUser("Region crop #" + crop350Num,
			"On CONFOCAL_MAX: drag the yellow " + side350Px + "×" + side350Px + " px square.\n" +
			"Target labels and ovals are visible on this image.\n" +
			"Click OK when the box is positioned.");

		selectImage(confocalMaxID);
		if (selectionType() != 0 && roiManager("count") > 0) {
			roiManager("Update");
		}
		boundsStr = getCropBoundsFromImage();
		if (boundsStr == "") boundsStr = getCropBoundsFromManager(0);
		if (boundsStr == "") {
			showMessage("Could not read the crop square (ROI count=" + roiManager("count") + ").\n" +
				"Drag the yellow rectangle in the ROI Manager list, then try again.");
			continue;
		}

		parts = split(boundsStr, ",");
		cropX = parseInt(parts[0]);
		cropY = parseInt(parts[1]);
		cropW = parseInt(parts[2]);
		cropH = parseInt(parts[3]);

		if (!getBoolean("Save region crop #" + crop350Num + " at (" + cropX + ", " + cropY + ")?")) {
			if (getBoolean("Finished with region crops?")) break;
			continue;
		}

		confocalCropPath = outputDir + baseName + "_intermediate_crop_" + pad2(crop350Num) + "_confocal.png";
		reflectionCropPath = outputDir + baseName + "_intermediate_crop_" + pad2(crop350Num) + "_reflection.png";
		duplicateCropToFile(confocalMaxID, cropX, cropY, cropW, cropH, confocalCropPath);
		duplicateCropToFile(reflectionMaxID, cropX, cropY, cropW, cropH, reflectionCropPath);
		appendRegionCropCsv(cropCsvPath, crop350Num, cropX, cropY, cropW, cropH,
			confocalCropPath, reflectionCropPath, cropCsvReady);
		cropCsvReady = 1;

		crop350Num++;
		if (!getBoolean("Add another region crop?")) break;
	}

	// ---- Scale calibration for target zoom scale bars (user-defined only) ----
	showMessage(
		"Target zoom exports include a 50 µm scale bar (bottom right).\n\n" +
		"Enter how many pixels match a known distance in µm.\n" +
		"Example: if 100 px = 10 µm, enter 100 and 10."
	);
	calPx = getNumber("Calibration — distance (pixels)", 100);
	calUm = getNumber("Calibration — known distance (µm)", 10);
	if (calPx <= 0 || calUm <= 0) exit("Calibration values must be greater than zero.");
	umPerPx = calUm / calPx;
	scaleBarPx = 50 / umPerPx;
	if (scaleBarPx > side80Px - 20) {
		exit("50 µm scale bar (" + d2s(scaleBarPx, 1) + " px) is wider than the " +
			side80Px + " px zoom crop.\nUse a finer calibration (more µm per pixel span).");
	}

	// ---- 4. Interactive 350 px target zooms on the spotted Z-slice (no marks) ----
	setOption("Show All", false);

	for (tIdx = 0; tIdx < targetCount; tIdx++) {
		tgtNum = tIdx + 1;
		sliceZ = targetZ[tIdx];

		roiManager("reset");
		selectImage(stackID);
		Overlay.remove;
		run("Select None");

		sliceViewID = buildConfocalSliceView(stackID, sliceZ, hasGreen);
		selectImage(sliceViewID);
		setTool("rectangle");
		centerSquareOnTarget(tIdx, side80Px);
		roiManager("Add");
		roiManager("Select", 0);

		if (roiManager("count") == 0) {
			showMessage("Could not create zoom rectangle. Try again.");
			close(sliceViewID);
			tIdx--;
			continue;
		}

		msg = "Target " + tgtNum + " of " + targetCount + "  (Z " + sliceZ + "/" + numZ + ")\n\n" +
			"Confocal view at the Z slice where this target was marked.\n" +
			"Drag the 350×350 px square if needed, then click OK.";
		waitForUser("Target zoom — target " + tgtNum, msg);

		selectImage(sliceViewID);
		if (selectionType() != 0 && roiManager("count") > 0) {
			roiManager("Update");
		}
		boundsStr = getCropBoundsFromImage();
		if (boundsStr == "") boundsStr = getCropBoundsFromManager(0);
		if (boundsStr == "") {
			showMessage("Could not read the zoom square. Try again.");
			close(sliceViewID);
			tIdx--;
			continue;
		}
		parts = split(boundsStr, ",");
		cropX = parseInt(parts[0]);
		cropY = parseInt(parts[1]);
		cropW = parseInt(parts[2]);
		cropH = parseInt(parts[3]);

		if (getBoolean("Save 350 px zoom PNG for target " + tgtNum + "?")) {
			zoomPath = outputDir + baseName + "_zoomed_target" + pad2(tgtNum) + "_z" + pad2(sliceZ) + "_confocal.png";
			saveZoomCropWithScaleBar(sliceViewID, cropX, cropY, cropW, cropH, zoomPath, calPx, calUm);
		}
		close(sliceViewID);
	}

	close(confocalMaxID);
	close(reflectionMaxID);
	selectImage(stackID);
	close();

	showMessage(
		"Done!\n\nSaved to:\n" + outputDir + "\n\n" +
		"- " + baseName + "_max_confocal.png\n" +
		"- " + baseName + "_max_reflection.png\n" +
		"- " + baseName + "_intermediate_crop_XX_confocal.png / _reflection.png\n" +
		"- " + baseName + "_intermediate_crops.csv\n" +
		"- " + baseName + "_zoomed_targetXX_zZZ_confocal.png (50 µm scale bar, saved targets only)"
	);
}

// ---- ROI overlay (reads globals set in main; nRois passed explicitly) ----

function burnStoredRoisOnImage(windowName, nRois) {
	Overlay.remove;
	setFont("SansSerif", 14, "bold");
	markCount = 0;

	for (i = 0; i < nRois; i++) {
		if (roiIsLabel[i] == 1) {
			makeText(roiText[i], roiX[i], roiY[i]);
			Overlay.addSelection("white", 1);
			markCount++;
		} else if (roiIsTarget[i] == 1) {
			makeOval(roiX[i], roiY[i], roiW[i], roiH[i]);
			Overlay.addSelection("yellow", 2);
			markCount++;
		}
	}

	run("Select None");
	if (markCount > 0) Overlay.flatten;
	rename(windowName);
	return getImageID();
}

function getTextRoiText() {
	txt = eval("js",
		"importClass(Packages.ij.gui.TextRoi);" +
		"var imp = IJ.getImage();" +
		"var roi = imp.getRoi();" +
		"if (roi != null && roi instanceof TextRoi) roi.getText(); else '';");
	if (txt == "") txt = Roi.getName();
	return txt;
}

function clearAllRoisFromImage(imageID) {
	selectImage(imageID);
	n = roiManager("count");
	for (i = n - 1; i >= 0; i--) {
		roiManager("Select", i);
		roiManager("Delete");
	}
	roiManager("reset");
	Overlay.remove;
	run("Select None");
}

function buildReflectionMaxProjection(sourceID, ch, zCount) {
	reflMaxID = maxProjectOneChannel(sourceID, ch, zCount);
	selectImage(reflMaxID);
	rename("REFLECTION_MAX");
	applyGrayscaleDisplay();
	return getImageID();
}

function applyGrayscaleDisplay() {
	if (is("composite")) run("Make Composite");
	run("8-bit");
	run("Grays");
}

function applyBlueDisplay() {
	if (is("composite")) run("Make Composite");
	run("8-bit");
	run("Blue");
	run("RGB Color");
}

function parseTargetDrawnZ(name) {
	marker = "_drawnZ";
	pos = indexOf(name, marker);
	if (pos < 0) return 1;
	return parseInt(substring(name, pos + 7, lengthOf(name)));
}

// ---- Max projections ----
// 2-channel stacks (blue + reflection): max-project each channel directly — no merge/composite.
// 3-channel stacks: merge blue + green for confocal; reflection is channel 3 alone.

function buildConfocalMaxProjection(sourceID, includeGreen, zCount) {
	if (!includeGreen) {
		confocalMaxID = maxProjectOneChannel(sourceID, 1, zCount);
		selectImage(confocalMaxID);
		applyBlueDisplay();
		rename("CONFOCAL_MAX");
		return getImageID();
	}

	blueMaxID = maxProjectOneChannel(sourceID, 1, zCount);
	selectImage(blueMaxID);
	rename("BLUE_MAX");
	applyGrayscaleDisplay();

	greenMaxID = maxProjectOneChannel(sourceID, 2, zCount);
	selectImage(greenMaxID);
	rename("GREEN_MAX");
	applyGrayscaleDisplay();

	run("Merge Channels...", "c2=GREEN_MAX c3=BLUE_MAX create");
	close("BLUE_MAX", true);
	close("GREEN_MAX", true);
	run("Make Composite");
	rename("CONFOCAL_MAX");
	return getImageID();
}

function buildConfocalSliceView(sourceID, sliceZ, includeGreen) {
	if (!includeGreen) {
		sliceViewID = extractOneChannelSlice(sourceID, 1, sliceZ);
		selectImage(sliceViewID);
		applyBlueDisplay();
		return getImageID();
	}

	blueTag = "BLUE_Z" + sliceZ;
	greenTag = "GREEN_Z" + sliceZ;

	blueSliceID = extractOneChannelSlice(sourceID, 1, sliceZ);
	selectImage(blueSliceID);
	rename(blueTag);

	greenSliceID = extractOneChannelSlice(sourceID, 2, sliceZ);
	selectImage(greenSliceID);
	rename(greenTag);

	run("Merge Channels...", "c2=" + greenTag + " c3=" + blueTag + " create");
	close(blueTag, true);
	close(greenTag, true);
	run("Make Composite");
	return getImageID();
}

function maxProjectOneChannel(sourceID, ch, zCount) {
	selectImage(sourceID);
	stackTag = "CH" + ch + "_STACK";
	projTag = "CH" + ch + "_MAX";

	eval("js",
		"importClass(Packages.ij.IJ);" +
		"importClass(Packages.ij.ImageStack);" +
		"importClass(Packages.ij.ImagePlus);" +
		"importClass(Packages.ij.plugin.ZProjector);" +
		"var src = IJ.getImage();" +
		"var w = src.getWidth();" +
		"var h = src.getHeight();" +
		"var ch = " + ch + ";" +
		"var zCount = " + zCount + ";" +
		"var stack = new ImageStack(w, h);" +
		"for (var z = 1; z <= zCount; z++) {" +
		"  src.setPosition(ch, z, 1);" +
		"  stack.addSlice(null, src.getProcessor().duplicate());" +
		"}" +
		"var stackImp = new ImagePlus('" + stackTag + "', stack);" +
		"var proj = ZProjector.run(stackImp, 'max');" +
		"if (proj != null) { proj.setTitle('" + projTag + "'); proj.show(); }");

	selectWindow(projTag);
	return getImageID();
}

function extractOneChannelSlice(sourceID, ch, sliceZ) {
	selectImage(sourceID);
	sliceTag = "CH" + ch + "_Z" + sliceZ;

	eval("js",
		"importClass(Packages.ij.IJ);" +
		"importClass(Packages.ij.ImagePlus);" +
		"var src = IJ.getImage();" +
		"var ch = " + ch + ";" +
		"var z = " + sliceZ + ";" +
		"src.setPosition(ch, z, 1);" +
		"var imp = new ImagePlus('" + sliceTag + "', src.getProcessor().duplicate());" +
		"imp.show();");

	selectWindow(sliceTag);
	return getImageID();
}

// ---- Cropping ----

function getCropBoundsFromImage() {
	return eval("js",
		"importClass(Packages.ij.IJ);" +
		"var imp = IJ.getImage();" +
		"(imp == null) ? '' : (function() {" +
		"  var roi = imp.getRoi();" +
		"  if (roi == null) return '';" +
		"  var b = roi.getBounds();" +
		"  return b.x + ',' + b.y + ',' + b.width + ',' + b.height;" +
		"})();");
}

function getCropBoundsFromManager(roiIdx) {
	return eval("js",
		"importClass(Packages.ij.plugin.frame.RoiManager);" +
		"var rm = RoiManager.getInstance();" +
		"(rm == null || rm.getCount() <= " + roiIdx + ") ? '' : " +
		"(function() {" +
		"  var roi = rm.getRoi(" + roiIdx + ");" +
		"  if (roi == null) return '';" +
		"  var b = roi.getBounds();" +
		"  return b.x + ',' + b.y + ',' + b.width + ',' + b.height;" +
		"})();");
}

function duplicateCropToFile(sourceID, bx, by, bw, bh, savePath) {
	selectImage(sourceID);
	makeRectangle(bx, by, bw, bh);
	run("Duplicate...", "duplicate");
	saveAs("PNG", savePath);
	close();
	selectImage(sourceID);
	run("Select None");
}

function saveZoomCropWithScaleBar(sourceID, bx, by, bw, bh, savePath, calPx, calUm) {
	selectImage(sourceID);
	makeRectangle(bx, by, bw, bh);
	run("Duplicate...", "duplicate");

	run("Set Scale...", "distance=" + d2s(calPx, 9) + " known=" + d2s(calUm, 9) + " unit=µm");
	run("Scale Bar...", "width=50 height=10 font=16 color=White location=[Lower Right] bold overlay");
	Overlay.flatten;
	saveAs("PNG", savePath);
	close();
	selectImage(sourceID);
	run("Select None");
}

function appendRegionCropCsv(csvPath, cropNum, bx, by, bw, bh, confocalPath, reflectionPath, headerDone) {
	csvHandle = File.open(csvPath);
	if (headerDone == 0) {
		print(csvHandle, "Crop,TopLeftX_px,TopLeftY_px,Width_px,Height_px,ConfocalPNG,ReflectionPNG");
	}
	print(csvHandle, "" + cropNum + "," + bx + "," + by + "," + bw + "," + bh + "," +
		File.getName(confocalPath) + "," + File.getName(reflectionPath));
	File.close(csvHandle);
}

function positionSquareRoi(sidePx) {
	imgW = getWidth();
	imgH = getHeight();
	if (sidePx > imgW || sidePx > imgH) {
		exit("Image is too small for a " + sidePx + " px square crop.\nImage size: " + imgW + " x " + imgH);
	}
	startX = round((imgW - sidePx) / 2);
	startY = round((imgH - sidePx) / 2);
	if (startX < 0) startX = 0;
	if (startY < 0) startY = 0;
	if (startX + sidePx > imgW) startX = imgW - sidePx;
	if (startY + sidePx > imgH) startY = imgH - sidePx;

	setTool("rectangle");
	makeRectangle(startX, startY, sidePx, sidePx);
}

function centerSquareOnTarget(tIdx, sidePx) {
	cx = targetCenterX[tIdx];
	cy = targetCenterY[tIdx];
	half = sidePx / 2;
	rx = round(cx - half);
	ry = round(cy - half);
	imgW = getWidth();
	imgH = getHeight();
	if (sidePx > imgW || sidePx > imgH) {
		exit("Image is too small for a " + sidePx + " px square crop.\nImage size: " + imgW + " x " + imgH);
	}
	if (rx < 0) rx = 0;
	if (ry < 0) ry = 0;
	if (rx + sidePx > imgW) rx = imgW - sidePx;
	if (ry + sidePx > imgH) ry = imgH - sidePx;

	setTool("rectangle");
	makeRectangle(rx, ry, sidePx, sidePx);
}

// ---- File helpers ----

function findCziFile(dir) {
	list = getFileList(dir);
	if (lengthOf(list) == 0) return "";

	foundPath = "";
	cziCount = 0;
	for (i = 0; i < lengthOf(list); i++) {
		if (endsWith(list[i], ".czi") || endsWith(list[i], ".CZI")) {
			cziCount++;
			foundPath = dir + list[i];
		}
	}

	if (cziCount == 1) return foundPath;

	for (i = 0; i < lengthOf(list); i++) {
		if (!endsWith(list[i], ".czi") && !endsWith(list[i], ".CZI")) continue;
		candidate = dir + list[i];
		candidateBase = getBaseName(candidate);
		zipCandidate = dir + candidateBase + "_targets.zip";
		if (File.exists(zipCandidate)) return candidate;
	}
	return foundPath;
}

function findTargetsZip(dir, baseName) {
	candidate = dir + baseName + "_targets.zip";
	if (File.exists(candidate)) return candidate;
	return "";
}

function getBaseName(path) {
	name = File.getName(path);
	dot = lastIndexOf(name, ".");
	if (dot > 0) name = substring(name, 0, dot);
	return name;
}

function pad2(n) {
	if (n < 10) return "0" + n;
	return "" + n;
}
