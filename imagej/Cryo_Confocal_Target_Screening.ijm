// Cryo-Confocal Target Screening Macro
// Fiji / ImageJ 1.54p (requires Bio-Formats)
//
// Workflow:
//   1. Choose output folder and import a .czi Z-stack (multi-channel supported)
//   2. Scroll Z, zoom, draw ovals around targets until finished
//   3. Saves max-projection PNG, ROI set, CSV, and a copy of the CZI

macro "Cryo Confocal Target Screening" {
	requires("1.52n");

	outputDir = getDirectory("Choose output folder");
	if (outputDir == "") exit("No output folder selected.");
	if (!endsWith(outputDir, File.separator)) outputDir = outputDir + File.separator;

	inputPath = File.openDialog("Select CZI file");
	if (inputPath == "") exit("No file selected.");

	baseName = getBaseName(inputPath);

	// Import as hyperstack; composite view for blue/green channels
	run("Bio-Formats Importer", "open=[" + inputPath + "] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT");

	stackID = getImageID();
	stackTitle = getTitle();
	getDimensions(impW, impH, numChannels, numZ, numFrames);
	totalZ = numZ;

	if (numChannels > 1) run("Make Composite");

	// Target storage (expandable up to maxTargets)
	maxTargets = 500;
	targetX = newArray(maxTargets);
	targetY = newArray(maxTargets);
	targetW = newArray(maxTargets);
	targetH = newArray(maxTargets);
	targetZ = newArray(maxTargets);
	nTargets = 0;

	roiManager("reset");
	roiManager("Associate", "true");
	setOption("Show All", true);
	setTool("oval");

	showMessage(
		"Target screening started.\n\n" +
		"Use the Z slider and zoom tools to find structures.\n" +
		"Draw ovals with the ellipse tool.\n\n" +
		"OK with an oval selected  = add target\n" +
		"OK without a selection    = finish and save\n\n" +
		"Do not click Cancel on the prompt (that stops the macro)."
	);

	// ---- Interactive target marking ----
	while (true) {
		updateScreeningDisplay();

		msg = "Targets marked: " + nTargets + "\n\n" +
			"1. Scroll to the correct Z slice\n" +
			"2. Zoom as needed\n" +
			"3. Draw an oval around the target\n\n" +
			"OK with oval selected  -> add target\n" +
			"OK without selection   -> finish and save";
		waitForUser("Add Target", msg);

		if (selectionType() == -1) {
			if (nTargets == 0) exit("No targets were marked.");
			if (getBoolean("Save outputs for " + nTargets + " target(s)?")) break;
			continue;
		}

		if (selectionType() != 1) {
			showMessage("Please draw an oval (ellipse) selection.\nOther selection types are not supported.");
			run("Select None");
			continue;
		}

		if (nTargets >= maxTargets) exit("Maximum number of targets (" + maxTargets + ") reached.");

		getSelectionBounds(x, y, w, h);
		Stack.getPosition(ch, zDrawn, f);

		targetX[nTargets] = x;
		targetY[nTargets] = y;
		targetW[nTargets] = w;
		targetH[nTargets] = h;
		targetZ[nTargets] = zDrawn;
		nTargets++;

		// Save oval ROI for collaborators (visible on all Z planes; z=0 means every slice)
		selectImage(stackID);
		makeOval(x, y, w, h);
		Roi.setPosition(1, 0, 1);
		Roi.setName("Target_" + pad2(nTargets) + "_drawnZ" + zDrawn);
		Roi.setStrokeColor(255, 255, 0);
		Roi.setStrokeWidth(2);
		roiManager("Add");

		// Label ROI so the user can see what is already marked while screening
		markText = buildTargetLabel(nTargets, zDrawn, totalZ);
		textX = x + w / 2 - getLabelOffset(nTargets, zDrawn, totalZ);
		textY = y - 6;
		if (textY < 2) textY = y + h + 14;
		setFont("SansSerif", 14, "bold");
		makeText(markText, textX, textY);
		Roi.setPosition(1, 0, 1);
		Roi.setName("Label_" + pad2(nTargets));
		Roi.setStrokeColor(255, 255, 255);
		roiManager("Add");

		run("Select None");
		updateScreeningDisplay();
	}

	// ---- Save outputs ----
	setBatchMode(true);
	selectImage(stackID);

	// Max-intensity projection with target marks (yellow ovals, white labels)
	projTitle = createMaxIntensityProjection(stackID, numChannels, totalZ);
	selectWindow(projTitle);
	prepareForAnnotation(projTitle, numChannels);
	drawTargetsOnImage();
	saveAs("PNG", outputDir + baseName + "_max_projection.png");
	close();

	// ROI set
	selectImage(stackID);
	roiManager("Save", outputDir + baseName + "_targets.zip");

	// CSV with pixel coordinates (written directly; avoids reserved names like "line")
	csvPath = outputDir + baseName + "_targets.csv";
	csvHandle = File.open(csvPath);
	print(csvHandle, "Target,Label,CenterX_px,CenterY_px,TopLeftX_px,TopLeftY_px,Width_px,Height_px,Z_drawn,Total_Z,ImageWidth_px,ImageHeight_px");
	for (tIdx = 0; tIdx < nTargets; tIdx++) {
		tgtNum = tIdx + 1;
		markText = buildTargetLabel(tgtNum, targetZ[tIdx], totalZ);
		centerX = targetX[tIdx] + targetW[tIdx] / 2;
		centerY = targetY[tIdx] + targetH[tIdx] / 2;
		print(csvHandle, "" + tgtNum + "," + markText + "," + d2s(centerX, 2) + "," + d2s(centerY, 2) + "," +
			d2s(targetX[tIdx], 0) + "," + d2s(targetY[tIdx], 0) + "," +
			d2s(targetW[tIdx], 0) + "," + d2s(targetH[tIdx], 0) + "," +
			d2s(targetZ[tIdx], 0) + "," + d2s(totalZ, 0) + "," + d2s(impW, 0) + "," + d2s(impH, 0));
	}
	File.close(csvHandle);

	// Copy original CZI so collaborators can reopen it with the ROI set
	File.copy(inputPath, outputDir + baseName + ".czi");

	setBatchMode(false);
	selectImage(stackID);
	updateScreeningDisplay();

	showMessage(
		"Done!\n\n" +
		"Saved to:\n" + outputDir + "\n\n" +
		"- " + baseName + ".czi (copy of input)\n" +
		"- " + baseName + "_targets.zip\n" +
		"- " + baseName + "_targets.csv\n" +
		"- " + baseName + "_max_projection.png"
	);
}

// Keep marked targets visible on every Z slice while screening
function updateScreeningDisplay() {
	selectImage(stackID);
	setOption("Show All", true);
	roiManager("Show All");
}

function createMaxIntensityProjection(sourceID, chCount, zCount) {
	// Max-project confocal channels only (exclude reflection when 3 channels)
	if (chCount >= 3) {
		blueMaxID = maxProjectOneChannel(sourceID, 1, zCount);
		selectImage(blueMaxID);
		rename("BLUE_MAX");
		greenMaxID = maxProjectOneChannel(sourceID, 2, zCount);
		selectImage(greenMaxID);
		rename("GREEN_MAX");
		run("Merge Channels...", "c2=GREEN_MAX c3=BLUE_MAX create");
		close("BLUE_MAX");
		close("GREEN_MAX");
		run("Make Composite");
		rename("MAXPROJ");
		return getTitle();
	}

	if (chCount == 2) {
		// 2 channels = blue (C1) + reflection (C2); confocal max is blue channel only
		maxProjectOneChannel(sourceID, 1, zCount);
		selectImage("CH1_MAX");
		run("8-bit");
		run("Blue");
		run("RGB Color");
		rename("MAXPROJ");
		return getTitle();
	}

	maxProjectOneChannel(sourceID, 1, zCount);
	rename("MAXPROJ");
	return getTitle();
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

function prepareForAnnotation(title, nChannels) {
	selectWindow(title);
	if (nChannels > 1) {
		if (!is("composite")) run("Make Composite");
	}
}

function drawTargetsOnImage() {
	Overlay.remove;
	setFont("SansSerif", 14, "bold");

	for (i = 0; i < nTargets; i++) {
		markText = buildTargetLabel(i + 1, targetZ[i], totalZ);
		textX = targetX[i] + targetW[i] / 2 - getLabelOffset(i + 1, targetZ[i], totalZ);
		textY = targetY[i] - 6;
		if (textY < 2) textY = targetY[i] + targetH[i] + 14;

		makeOval(targetX[i], targetY[i], targetW[i], targetH[i]);
		Overlay.addSelection("yellow", 2);

		makeText(markText, textX, textY);
		Overlay.addSelection("white", 1);
	}

	run("Select None");
	Overlay.flatten;
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

function buildTargetLabel(targetNum, zDrawn, zTotal) {
	// ImageJ parses bare "/" in "a + "/" + b" as division, not a slash character
	zSep = "/";
	return "" + targetNum + " (z: " + zDrawn + zSep + zTotal + ")";
}

function getLabelOffset(targetNum, zDrawn, zTotal) {
	txt = buildTargetLabel(targetNum, zDrawn, zTotal);
	// Rough centering for overlay/saved text (approx. 7 px per character at 14 pt)
	return (lengthOf(txt) * 7) / 2;
}
