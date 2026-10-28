// follicle_only_measurement.ijm
//
// Fiji / ImageJ macro language
//
// Simplified follicle-only workflow derived from
// multi_oocyte_largest_slice_measurement.ijm
//
// Measures follicle area, major/minor axis, and eccentricity on the
// clearest z-slice. Records top/bottom z with visibility status and
// computes z-extent only when both ends are visible.
//
// Display prep: per-channel auto/manual brightness-contrast, then user
// selects which channels to merge (any combination).
//
// Input: experiment folder containing one subfolder per z-stack (TIFFs +
// MetaData per subfolder). All stacks append to one CSV. You choose which
// stack to open first and which to open next after each stack is finished.

requires("1.53");

// --------------------------------------------------
// Settings
// --------------------------------------------------
sep = ";";
numDigits = 6;

lastCalMode = "direct_pixel_size";
lastDirectPxW = 1.0;
lastDirectPxH = 1.0;
lastFullWidth = 512.0;
lastFullHeight = 512.0;
lastPixelsX = 512;
lastPixelsY = 512;
lastZStep = 2.0;
lastUnit = "um";
lastDoDisplayPrep = true;
lastUseAutoBC = true;

maxStackJobs = 500;
jobCount = 0;
jobLabel = newArray(maxStackJobs);
jobInputDir = newArray(maxStackJobs);
jobFileName = newArray(maxStackJobs);
jobFolderLabel = newArray(maxStackJobs);

// --------------------------------------------------
// 1) Select input/output
// --------------------------------------------------
experimentDir = getDirectory("Choose experiment folder containing z-stack subfolders");
if (experimentDir == "" || experimentDir == "null")
    exit("No experiment folder selected.");

if (!endsWith(experimentDir, File.separator))
    experimentDir = experimentDir + File.separator;

jobCount = 0;
folderPaths = collectStackFolderPaths(experimentDir);
for (sf = 0; sf < folderPaths.length; sf++) {
    dirPath = folderPaths[sf];
    if (!endsWith(dirPath, File.separator))
        dirPath = dirPath + File.separator;

    folderLabel = getFolderBaseName(dirPath);
    fileList = getFileList(dirPath);

    for (i = 0; i < fileList.length; i++) {
        name = fileList[i];
        if (isListEntryDirectory(name))
            continue;
        if (!endsWith(name, ".tif") && !endsWith(name, ".tiff"))
            continue;
        if (isLikelyCompanionSlice(name))
            continue;

        if (jobCount >= maxStackJobs)
            exit("Too many stack entry files (max " + maxStackJobs + ").");

        jobLabel[jobCount] = folderLabel + " / " + name;
        jobInputDir[jobCount] = dirPath;
        jobFileName[jobCount] = name;
        jobFolderLabel[jobCount] = folderLabel;
        jobCount++;
    }
}

if (jobCount < 1)
    exit("No stack entry files found in:\n" + experimentDir +
         "\n\nEach subfolder should contain *_z00.tif (and companion slices) and optionally MetaData\\.");

outputDir = getDirectory("Choose folder to save CSV and ROI files");
if (outputDir == "" || outputDir == "null")
    exit("No output folder selected.");

Dialog.create("Output settings");
Dialog.addString("CSV filename:", "follicle_only_measurements.csv");
Dialog.show();
csvName = Dialog.getString();
if (csvName == "") csvName = "follicle_only_measurements.csv";

outputCSV = outputDir + csvName;

roiDir = outputDir + "Saved_ROIs" + File.separator;
File.makeDirectory(roiDir);

header =
"StackSubfolder" + sep +
"Stack" + sep +
"FollicleID" + sep +
"FollicleBottomZ" + sep +
"FollicleTopZ" + sep +
"FollicleBottomVisibility" + sep +
"FollicleTopVisibility" + sep +
"FollicleBottomTopLabel" + sep +
"FollicleZExtent" + sep +
"FollicleZExtentUncertainty" + sep +
"FollicleROIZ" + sep +
"TotalZSlices" + sep +
"FollicleROIZLabel" + sep +
"PixelWidth" + sep +
"PixelHeight" + sep +
"ZStep" + sep +
"Unit" + sep +
"DisplayedChannels" + sep +
"CalibrationMode" + sep +
"MetadataSourceFile" + sep +
"MetadataLoadedAutomatically" + sep +
"MetadataPixelWidth" + sep +
"MetadataPixelHeight" + sep +
"MetadataZStep" + sep +
"MetadataUnit" + sep +
"MetadataPixelsX" + sep +
"MetadataPixelsY" + sep +
"MetadataZSlices" + sep +
"FullImageWidth" + sep +
"FullImageHeight" + sep +
"PixelsX" + sep +
"PixelsY" + sep +
"FollicleArea" + sep +
"FollicleMajor" + sep +
"FollicleMinor" + sep +
"FollicleEccentricity" + sep +
"FollicleType" + sep +
"GranulosaShape" + sep +
"PartialFollicle" + sep +
"WouldTarget" + sep +
"SignalQuality" + sep +
"BoundaryConfidence" + sep +
"ExcludeFromMainAnalysis" + sep +
"Notes" + sep +
"FollicleROIFile" + sep +
"ROIZipFile\n";

nStacksAlreadyInCsv = 0;

run("Set Measurements...", "area fit redirect=None decimal=6");

setBatchMode(false);

completedKeys = "|";
if (File.exists(outputCSV)) {
    csvText = File.openAsString(outputCSV);
    if (csvText != "") {
        csvLines = split(csvText, "\n");
        for (li = 1; li < csvLines.length; li++) {
            csvLine = csvLines[li];
            csvLine = replace(csvLine, "\r", "");
            if (csvLine == "")
                continue;

            csvParts = split(csvLine, ";");
            if (csvParts.length < 2)
                continue;

            csvSubfolder = trimCsvField(csvParts[0]);
            csvStackFile = trimCsvField(csvParts[1]);

            for (j = 0; j < jobCount; j++) {
                if (jobFolderLabel[j] == csvSubfolder && jobFileName[j] == csvStackFile) {
                    csvKey = jobFolderLabel[j] + "||" + jobFileName[j];
                    csvToken = "|" + csvKey + "|";
                    if (indexOf(completedKeys, csvToken) == -1) {
                        completedKeys = completedKeys + csvKey + "|";
                        nStacksAlreadyInCsv++;
                    }
                }
            }
        }
    }
} else {
    File.saveString(header, outputCSV);
}

if (countRemainingStackJobs() < 1)
    exit("All stack entry files are already recorded in:\n" + outputCSV);

showMessage("Stacks found",
    "Experiment: " + experimentDir + "\n\n" +
    "Stack entry files found: " + jobCount + "\n" +
    "Already in CSV (skipped): " + nStacksAlreadyInCsv + "\n" +
    "Still available: " + countRemainingStackJobs() + "\n\n" +
    "Choose which stack to open first.\n" +
    "After each stack, choose the next one or end the session.\n\n" +
    "All measurements append to:\n" + outputCSV);

jobIdx = pickStackJob("Select starting stack");
if (jobIdx < 0)
    exit("No stack selected.");

jobsDoneThisSession = 0;

// --------------------------------------------------
// 2) Process user-selected stacks until session ends
// --------------------------------------------------
while (jobIdx >= 0) {

    inputDir = jobInputDir[jobIdx];
    if (!endsWith(inputDir, File.separator))
        inputDir = inputDir + File.separator;

    name = jobFileName[jobIdx];
    stackFolderLabel = jobFolderLabel[jobIdx];
    stackRoiDir = roiDir + sanitizeName(stackFolderLabel) + File.separator;
    File.makeDirectory(stackRoiDir);

    fullPath = inputDir + name;
    stackBase = getSeriesBase(stripExtension(name));

    while (nImages > 0) {
        selectImage(nImages);
        close();
    }

    bfOptions =
        "open=[" + fullPath + "] " +
        "autoscale " +
        "group_files " +
        "color_mode=Composite " +
        "view=Hyperstack " +
        "stack_order=XYCZT";

    run("Bio-Formats Importer", bfOptions);

    if (nImages < 1)
        exit("Bio-Formats failed to open: " + name);

    rawTitle = stackBase + "_raw";
    run("Rename...", "title=[" + rawTitle + "]");

    getDimensions(imgW, imgH, nChannels, totalZSlices, nFrames);
    if (nChannels < 1) nChannels = 1;
    if (totalZSlices < 1) totalZSlices = 1;
    if (nFrames < 1) nFrames = 1;

    metadataSourceFile = inputDir + "MetaData" + File.separator + stackBase + "_Properties.xml";
    metadataLoadedAutomatically = "no";

    metadataPxW = -1;
    metadataPxH = -1;
    metadataZStep = -1;
    metadataUnit = "";
    metadataPixelsX = -1;
    metadataPixelsY = -1;
    metadataZSlices = -1;

    if (File.exists(metadataSourceFile)) {
        xmlText = File.openAsString(metadataSourceFile);

        if (xmlText != "" && indexOf(xmlText, "<Dimensions") != -1) {
            xTag = getDimensionTag(xmlText, "X");
            yTag = getDimensionTag(xmlText, "Y");
            zTag = getDimensionTag(xmlText, "Z");

            if (xTag != "" && yTag != "" && zTag != "") {
                metadataPxW = parseNumberSmart(getXmlAttr(xTag, "Voxel"));
                metadataPxH = parseNumberSmart(getXmlAttr(yTag, "Voxel"));
                metadataZStep = parseNumberSmart(getXmlAttr(zTag, "Voxel"));

                metadataPixelsX = floor(parseNumberSmart(getXmlAttr(xTag, "NumberOfElements")));
                metadataPixelsY = floor(parseNumberSmart(getXmlAttr(yTag, "NumberOfElements")));
                metadataZSlices = floor(parseNumberSmart(getXmlAttr(zTag, "NumberOfElements")));

                metadataUnit = normalizeUnit(getXmlAttr(xTag, "Unit"));

                if (metadataPxW > 0 && metadataPxH > 0 && metadataZStep > 0 &&
                    metadataPixelsX > 0 && metadataPixelsY > 0 && metadataZSlices > 0 &&
                    metadataUnit != "") {
                    metadataLoadedAutomatically = "yes";
                }
            }
        }
    }

    if (metadataLoadedAutomatically == "yes") {
        pxW = metadataPxW;
        pxH = metadataPxH;
        zStep = metadataZStep;
        unit = metadataUnit;

        if (unit == "") unit = "um";

        fullWidthUsed = pxW * metadataPixelsX;
        fullHeightUsed = pxH * metadataPixelsY;
        pixelsXUsed = metadataPixelsX;
        pixelsYUsed = metadataPixelsY;

        calMode = "metadata_auto";
        doDisplayPrep = lastDoDisplayPrep;

        lastDirectPxW = pxW;
        lastDirectPxH = pxH;
        lastFullWidth = fullWidthUsed;
        lastFullHeight = fullHeightUsed;
        lastPixelsX = pixelsXUsed;
        lastPixelsY = pixelsYUsed;
        lastZStep = zStep;
        lastUnit = unit;

    } else {
        getVoxelSize(importPxW, importPxH, importZStep, importUnit);
        if (importUnit == "") importUnit = "unknown";
        if (importPxW <= 0) importPxW = 1;
        if (importPxH <= 0) importPxH = 1;
        if (importZStep <= 0) importZStep = 1;

        if (jobsDoneThisSession == 0) {
            if (importPxW > 0 && importPxW < 100) lastDirectPxW = importPxW;
            if (importPxH > 0 && importPxH < 100) lastDirectPxH = importPxH;
            if (importZStep > 0 && importZStep < 1000) lastZStep = importZStep;
            if (importUnit != "unknown" && importUnit != "inch" && importUnit != "inches") lastUnit = importUnit;
        }

        if (lastPixelsX <= 0) lastPixelsX = imgW;
        if (lastPixelsY <= 0) lastPixelsY = imgH;
        if (lastFullWidth <= 0) lastFullWidth = lastDirectPxW * imgW;
        if (lastFullHeight <= 0) lastFullHeight = lastDirectPxH * imgH;

        Dialog.create("Calibration fallback for " + name);
        Dialog.addMessage("Automatic metadata loading failed.\nUsing manual calibration dialog.");
        Dialog.addMessage("Expected metadata file:\n" + metadataSourceFile);
        modeOptions = newArray("direct_pixel_size", "full_image_size");
        Dialog.addChoice("Calibration entry mode:", modeOptions, lastCalMode);

        Dialog.addNumber("Direct pixel width:", lastDirectPxW);
        Dialog.addNumber("Direct pixel height:", lastDirectPxH);

        Dialog.addNumber("Full image width:", lastFullWidth);
        Dialog.addNumber("Full image height:", lastFullHeight);
        Dialog.addNumber("Pixels in X:", imgW);
        Dialog.addNumber("Pixels in Y:", imgH);

        Dialog.addNumber("Z-step / voxel depth:", lastZStep);
        Dialog.addString("Unit:", lastUnit);
        Dialog.addCheckbox("Prepare display (split channels, adjust red/blue, remove green)", lastDoDisplayPrep);
        Dialog.show();

        calMode = Dialog.getChoice();

        directPxWInput = Dialog.getNumber();
        directPxHInput = Dialog.getNumber();

        fullWidthInput = Dialog.getNumber();
        fullHeightInput = Dialog.getNumber();
        pixelsXInput = Dialog.getNumber();
        pixelsYInput = Dialog.getNumber();

        zStep = Dialog.getNumber();
        unit = Dialog.getString();
        doDisplayPrep = Dialog.getCheckbox();

        if (unit == "") unit = "um";

        if (zStep <= 0)
            exit("Z-step must be > 0.");

        if (calMode == "direct_pixel_size") {
            pxW = directPxWInput;
            pxH = directPxHInput;
            fullWidthUsed = pxW * imgW;
            fullHeightUsed = pxH * imgH;
            pixelsXUsed = imgW;
            pixelsYUsed = imgH;

            if (pxW <= 0 || pxH <= 0)
                exit("Direct pixel width and height must be > 0.");
        } else {
            fullWidthUsed = fullWidthInput;
            fullHeightUsed = fullHeightInput;
            pixelsXUsed = pixelsXInput;
            pixelsYUsed = pixelsYInput;

            if (fullWidthUsed <= 0 || fullHeightUsed <= 0 || pixelsXUsed <= 0 || pixelsYUsed <= 0)
                exit("Full image width/height and pixels in X/Y must be > 0.");

            pxW = fullWidthUsed / pixelsXUsed;
            pxH = fullHeightUsed / pixelsYUsed;

            if (pxW <= 0 || pxH <= 0)
                exit("Computed pixel size is invalid.");
        }

        lastCalMode = calMode;
        lastDirectPxW = pxW;
        lastDirectPxH = pxH;
        lastFullWidth = fullWidthUsed;
        lastFullHeight = fullHeightUsed;
        lastPixelsX = pixelsXUsed;
        lastPixelsY = pixelsYUsed;
        lastZStep = zStep;
        lastUnit = unit;
        lastDoDisplayPrep = doDisplayPrep;
    }

    setVoxelSize(pxW, pxH, zStep, unit);

    displayedChannels = "all";

    if (doDisplayPrep) {
        displayedChannels = prepareDisplayStack(rawTitle, stackBase, name, nChannels);

        browseMsg = "Browse the prepared stack and count follicles.\nClick OK when ready to continue.";
        if (nChannels > 1)
            browseMsg = "Selected channels were adjusted and merged.\n" + browseMsg;
        waitForUser("Image ready for " + name, browseMsg);

    } else {
        selectWindow(rawTitle);
        run("Rename...", "title=[" + stackBase + "]");
        if (nChannels == 3)
            displayedChannels = "111";
        else if (nChannels == 2)
            displayedChannels = "11";
        else
            displayedChannels = "all";
    }

    getDimensions(imgW, imgH, nChannels, totalZSlices, nFrames);
    if (nChannels < 1) nChannels = 1;
    if (totalZSlices < 1) totalZSlices = 1;
    if (nFrames < 1) nFrames = 1;

    run("Show Overlay");

    if (totalZSlices == 1) {
        waitForUser("Check stack import",
            "File: " + name + "\n\n" +
            "The current image reports only 1 Z slice.\n" +
            "If this should be a z-stack, the files may not be grouped correctly,\n" +
            "or this series may truly be only one plane.\n\n" +
            "Click OK to continue anyway.");
    }

    Dialog.create("Follicles in stack");
    Dialog.addMessage("Stack: " + name);
    Dialog.addNumber("Number of follicles to measure in this stack:", 1);
    Dialog.show();
    nFollicles = floor(Dialog.getNumber());

    if (nFollicles < 1) {
        waitForUser("Review image before continuing",
            "Stack: " + name + "\n\n" +
            "The image remains open.\n" +
            "Take any screenshot you want.\n" +
            "Then manually close the image window and click OK to continue.");
        while (isOpen(stackBase))
            waitForUser("Waiting for manual close",
                "Please manually close the image window titled:\n" + stackBase + "\n\nThen click OK.");
    } else {

    // ----------------------------------------------
    // 3) Measure each follicle in current stack
    // ----------------------------------------------
    for (ff = 1; ff <= nFollicles; ff++) {

        Dialog.create("Follicle ID");
        Dialog.addMessage("Stack: " + name);
        Dialog.addString("Enter ID for follicle " + ff + ":", "F" + ff);
        Dialog.show();
        follicleID = Dialog.getString();
        follicleID = sanitizeName(follicleID);
        if (follicleID == "") follicleID = "F" + ff;

        roiManager("reset");
        run("Clear Results");

        showMessage("Mark follicle extent for " + follicleID,
            "Stack: " + name + "\n" +
            "Follicle ID: " + follicleID + "\n\n" +
            "First, mark the BOTTOM and TOP z-planes of the follicle.");

        waitForUser("Select BOTTOM z-plane for " + follicleID,
            "Scroll to the BOTTOM of the follicle.\n" +
            "Click OK when you are on the best-estimate bottom plane.");
        follicleBottomZ_raw = getCurrentZSlice();

        visibilityOptions = newArray("visible", "out_of_focus");
        Dialog.create("Bottom visibility for " + follicleID);
        Dialog.addMessage("Bottom z-slice: " + follicleBottomZ_raw + "/" + totalZSlices);
        Dialog.addChoice("Is the bottom of the follicle visible on this plane?", visibilityOptions, "visible");
        Dialog.show();
        follicleBottomVisibility = Dialog.getChoice();

        waitForUser("Select TOP z-plane for " + follicleID,
            "Scroll to the TOP of the follicle.\n" +
            "Click OK when you are on the best-estimate top plane.");
        follicleTopZ_raw = getCurrentZSlice();

        Dialog.create("Top visibility for " + follicleID);
        Dialog.addMessage("Top z-slice: " + follicleTopZ_raw + "/" + totalZSlices);
        Dialog.addChoice("Is the top of the follicle visible on this plane?", visibilityOptions, "visible");
        Dialog.show();
        follicleTopVisibility = Dialog.getChoice();

        if (follicleBottomZ_raw <= follicleTopZ_raw) {
            follicleBottomZ = follicleBottomZ_raw;
            follicleTopZ = follicleTopZ_raw;
        } else {
            follicleBottomZ = follicleTopZ_raw;
            follicleTopZ = follicleBottomZ_raw;
        }

        follicleBottomTopLabel = "" + follicleBottomZ + "-" + follicleTopZ;

        if (follicleBottomVisibility == "visible" && follicleTopVisibility == "visible") {
            follicleZExtent = absVal(follicleTopZ - follicleBottomZ) * zStep;
            follicleZExtentUncertainty = zStep;
        } else {
            follicleZExtent = -1;
            follicleZExtentUncertainty = -1;
        }

        if (follicleBottomVisibility != "visible" || follicleTopVisibility != "visible")
            partialFollicle = "yes";
        else
            partialFollicle = "no";

        waitForUser("Choose FOLLICLE ROI z-plane for " + follicleID,
            "Scroll to the z-plane where the FOLLICLE cross-section is clearest.\n" +
            "Click OK when ready.");
        follicleROIZ = getCurrentZSlice();
        follicleROIZLabel = "" + follicleROIZ + "/" + totalZSlices;
        setCurrentZSlice(follicleROIZ);

        run("Select None");
        waitForUser("Draw FOLLICLE ROI: " + follicleID,
            "Stack: " + name + "\n" +
            "Follicle ID: " + follicleID + "\n" +
            "Follicle ROI z-slice: " + follicleROIZLabel + "\n\n" +
            "Draw the ROI around the ENTIRE FOLLICLE.\n" +
            "Then click OK.");

        if (selectionType() == -1) {
            showMessage("Skipped", "No follicle ROI drawn for " + follicleID + ". This follicle will be skipped.");
            run("Select None");
            roiManager("reset");
            continue;
        }

        getSelectionBounds(markerBX, markerBY, markerBW, markerBH);
        markerXPx = markerBX + markerBW / 2.0;
        markerYPx = markerBY + markerBH / 2.0;

        roiManager("Add");
        follicleIndex = roiManager("count") - 1;
        roiManager("Select", follicleIndex);
        roiManager("Rename", "FOLLICLE_" + stackBase + "_" + follicleID + "_z" + follicleROIZ + "_of_" + totalZSlices);

        run("Measure");
        follicleRow = nResults - 1;

        follicleArea  = getResult("Area", follicleRow);
        follicleMajor = getResult("Major", follicleRow);
        follicleMinor = getResult("Minor", follicleRow);

        if (follicleMajor > 0 && follicleMinor >= 0)
            follicleEcc = sqrt(1 - (follicleMinor * follicleMinor) / (follicleMajor * follicleMajor));
        else
            follicleEcc = NaN;

        Dialog.create("Follicle classification");
        Dialog.addMessage("Stack: " + name + "\nFollicle ID: " + follicleID +
                          "\nBottom z: " + follicleBottomZ + " (" + follicleBottomVisibility + ")" +
                          "\nTop z: " + follicleTopZ + " (" + follicleTopVisibility + ")" +
                          "\nFollicle ROI z: " + follicleROIZLabel +
                          "\nPartial follicle (auto): " + partialFollicle);

        classOptions = newArray("flat_primordial", "cuboidal_primordial", "primary", "unclear");
        granulosaOptions = newArray("flat", "cuboidal", "missing_few_flat", "missing_few_cuboidal", "multilayer");
        targetOptions = newArray("yes", "no");
        signalOptions = newArray("good", "bad");
        boundaryOptions = newArray("high", "medium", "low");
        excludeOptions = newArray("no", "yes");

        Dialog.addChoice("Follicle type:", classOptions, "unclear");
        Dialog.addChoice("Shape of granulosa cells:", granulosaOptions, "flat");
        Dialog.addChoice("Would you target:", targetOptions, "yes");
        Dialog.addChoice("Signal quality:", signalOptions, "good");
        Dialog.addChoice("Boundary confidence:", boundaryOptions, "high");
        Dialog.addChoice("Exclude from main analysis:", excludeOptions, "no");
        Dialog.addString("Notes:", "");
        Dialog.show();

        follicleType = Dialog.getChoice();
        granulosaShape = Dialog.getChoice();
        wouldTarget = Dialog.getChoice();
        signalQuality = Dialog.getChoice();
        boundaryConfidence = Dialog.getChoice();
        excludeFromMainAnalysis = Dialog.getChoice();
        notesText = Dialog.getString();

        follicleROIFile = stackRoiDir + stackBase + "_" + follicleID + "_z" + follicleROIZ + "_follicle.roi";
        roiZipFile      = stackRoiDir + stackBase + "_" + follicleID + "_ROIs.zip";

        roiManager("Select", follicleIndex);
        saveAs("Selection", follicleROIFile);

        roiManager("Deselect");
        roiManager("Save", roiZipFile);

        addPersistentMarkerAllZSlices(markerXPx, markerYPx, totalZSlices, follicleID);

        line = "";
        line = line + csvSafe(stackFolderLabel) + sep;
        line = line + csvSafe(name) + sep;
        line = line + csvSafe(follicleID) + sep;
        line = line + intToStr(follicleBottomZ) + sep;
        line = line + intToStr(follicleTopZ) + sep;
        line = line + csvSafe(follicleBottomVisibility) + sep;
        line = line + csvSafe(follicleTopVisibility) + sep;
        line = line + csvSafe(follicleBottomTopLabel) + sep;
        line = line + numToStr(follicleZExtent) + sep;
        line = line + numToStr(follicleZExtentUncertainty) + sep;
        line = line + intToStr(follicleROIZ) + sep;
        line = line + intToStr(totalZSlices) + sep;
        line = line + csvSafe(follicleROIZLabel) + sep;
        line = line + numToStr(pxW) + sep;
        line = line + numToStr(pxH) + sep;
        line = line + numToStr(zStep) + sep;
        line = line + csvSafe(unit) + sep;
        line = line + csvSafe(displayedChannels) + sep;
        line = line + csvSafe(calMode) + sep;
        line = line + csvSafe(metadataSourceFile) + sep;
        line = line + csvSafe(metadataLoadedAutomatically) + sep;
        line = line + numToStr(metadataPxW) + sep;
        line = line + numToStr(metadataPxH) + sep;
        line = line + numToStr(metadataZStep) + sep;
        line = line + csvSafe(metadataUnit) + sep;
        line = line + intToStr(metadataPixelsX) + sep;
        line = line + intToStr(metadataPixelsY) + sep;
        line = line + intToStr(metadataZSlices) + sep;
        line = line + numToStr(fullWidthUsed) + sep;
        line = line + numToStr(fullHeightUsed) + sep;
        line = line + intToStr(pixelsXUsed) + sep;
        line = line + intToStr(pixelsYUsed) + sep;
        line = line + numToStr(follicleArea) + sep;
        line = line + numToStr(follicleMajor) + sep;
        line = line + numToStr(follicleMinor) + sep;
        line = line + numToStr(follicleEcc) + sep;
        line = line + csvSafe(follicleType) + sep;
        line = line + csvSafe(granulosaShape) + sep;
        line = line + csvSafe(partialFollicle) + sep;
        line = line + csvSafe(wouldTarget) + sep;
        line = line + csvSafe(signalQuality) + sep;
        line = line + csvSafe(boundaryConfidence) + sep;
        line = line + csvSafe(excludeFromMainAnalysis) + sep;
        line = line + csvSafe(notesText) + sep;
        line = line + csvSafe(follicleROIFile) + sep;
        line = line + csvSafe(roiZipFile) + "\n";

        File.append(line, outputCSV);

        run("Clear Results");
        run("Select None");
        roiManager("reset");
        run("Show Overlay");
    }

    waitForUser("Take screenshot / review marked follicles",
        "Stack: " + name + "\n\n" +
        "The annotated image remains open so you can take a screenshot with all markers and names.\n\n" +
        "After you are done, manually close the image window titled:\n" + stackBase + "\n\n" +
        "Then click OK to continue.");

    while (isOpen(stackBase)) {
        waitForUser("Waiting for manual close",
            "Please manually close the image window titled:\n" + stackBase + "\n\n" +
            "Then click OK to continue to the next stack.");
    }

    } // end nFollicles > 0

    key = jobFolderLabel[jobIdx] + "||" + jobFileName[jobIdx];
    token = "|" + key + "|";
    if (indexOf(completedKeys, token) == -1)
        completedKeys = completedKeys + key + "|";

    jobsDoneThisSession++;

    if (countRemainingStackJobs() < 1)
        break;

    jobIdx = pickStackJob("Select next stack");

} // end stack session loop

showMessage("Done",
    "Measurement finished.\n\nResults saved to:\n" + outputCSV + "\n\nROI files saved in:\n" + roiDir);


// --------------------------------------------------
// Helper functions
// --------------------------------------------------

function countRemainingStackJobs() {
    n = 0;
    for (j = 0; j < jobCount; j++) {
        key = jobFolderLabel[j] + "||" + jobFileName[j];
        if (indexOf(completedKeys, "|" + key + "|") == -1)
            n = n + 1;
    }
    return n;
}

function pickStackJob(prompt) {
    nRemain = countRemainingStackJobs();
    if (nRemain < 1)
        return -1;

    choices = newArray(nRemain + 1);
    choiceJobIdx = newArray(nRemain + 1);
    c = 0;

    for (j = 0; j < jobCount; j++) {
        key = jobFolderLabel[j] + "||" + jobFileName[j];
        if (indexOf(completedKeys, "|" + key + "|") == -1) {
            choices[c] = jobLabel[j];
            choiceJobIdx[c] = j;
            c = c + 1;
        }
    }

    choices[c] = "--- End session ---";
    choiceJobIdx[c] = -1;

    Dialog.create(prompt);
    Dialog.addMessage("Remaining stacks: " + nRemain + " of " + jobCount);
    Dialog.addChoice("Stack to open:", choices, choices[0]);
    Dialog.show();
    selected = Dialog.getChoice();

    for (ci = 0; ci < choices.length; ci++) {
        if (choices[ci] == selected)
            return choiceJobIdx[ci];
    }
    return -1;
}

function trimCsvField(s) {
    return replace("" + s, "\r", "");
}

function collectStackFolderPaths(experimentDir) {
    if (!endsWith(experimentDir, File.separator))
        experimentDir = experimentDir + File.separator;

    list = getFileList(experimentDir);
    nFound = 0;

    for (i = 0; i < list.length; i++) {
        if (!isListEntryDirectory(list[i]))
            continue;
        entryName = stripListEntrySlash(list[i]);
        if (shouldSkipStackFolder(entryName))
            continue;
        subPath = experimentDir + list[i];
        if (hasProcessableStackEntry(subPath))
            nFound++;
    }

    if (nFound < 1) {
        if (hasProcessableStackEntry(experimentDir)) {
            paths = newArray(1);
            paths[0] = experimentDir;
            return paths;
        }
        return newArray(0);
    }

    paths = newArray(nFound);
    idx = 0;
    for (i = 0; i < list.length; i++) {
        if (!isListEntryDirectory(list[i]))
            continue;
        entryName = stripListEntrySlash(list[i]);
        if (shouldSkipStackFolder(entryName))
            continue;
        subPath = experimentDir + list[i];
        if (hasProcessableStackEntry(subPath)) {
            paths[idx] = subPath;
            idx++;
        }
    }
    return paths;
}

function hasProcessableStackEntry(dirPath) {
    if (!endsWith(dirPath, File.separator))
        dirPath = dirPath + File.separator;

    list = getFileList(dirPath);
    for (i = 0; i < list.length; i++) {
        if (isListEntryDirectory(list[i]))
            continue;
        name = list[i];
        if (!endsWith(name, ".tif") && !endsWith(name, ".tiff"))
            continue;
        if (isLikelyCompanionSlice(name))
            continue;
        return true;
    }
    return false;
}

function isListEntryDirectory(listEntry) {
    return endsWith(listEntry, "/") || endsWith(listEntry, "\\");
}

function stripListEntrySlash(listEntry) {
    n = lengthOf(listEntry);
    if (n < 1)
        return listEntry;
    last = substring(listEntry, n - 1, n);
    if (last == "/" || last == "\\")
        return substring(listEntry, 0, n - 1);
    return listEntry;
}

function shouldSkipStackFolder(folderName) {
    if (folderName == "MetaData")
        return true;
    if (folderName == "Saved_ROIs")
        return true;
    return false;
}

function getFolderBaseName(dirPath) {
    p = dirPath;
    while (lengthOf(p) > 0 && (endsWith(p, "/") || endsWith(p, "\\")))
        p = substring(p, 0, lengthOf(p) - 1);

    sepPos = lastIndexOf(p, File.separator);
    if (sepPos == -1)
        sepPos = lastIndexOf(p, "/");
    if (sepPos == -1)
        sepPos = lastIndexOf(p, "\\");
    if (sepPos == -1)
        return p;

    name = substring(p, sepPos + 1, lengthOf(p));
    while (lengthOf(name) > 0 && (endsWith(name, "/") || endsWith(name, "\\")))
        name = substring(name, 0, lengthOf(name) - 1);
    return name;
}

function prepareDisplayStack(rawTitle, stackBase, stackName, nChannels) {
    useAutoBC = lastUseAutoBC;
    Dialog.create("Display preparation for " + stackName);
    Dialog.addCheckbox("Auto brightness/contrast for each channel (Enhance Contrast)", lastUseAutoBC);
    Dialog.show();
    useAutoBC = Dialog.getCheckbox();
    lastUseAutoBC = useAutoBC;

    if (nChannels <= 1) {
        selectWindow(rawTitle);
        run("Rename...", "title=[" + stackBase + "]");
        adjustChannelStack(stackBase, stackName, "single channel", useAutoBC);
        return "all";
    }

    selectWindow(rawTitle);
    run("Split Channels");
    closeIfOpen(rawTitle);

    channelTitles = newArray(nChannels + 1);
    keepChannel = newArray(nChannels + 1);

    for (ch = 1; ch <= nChannels; ch++) {
        channelTitles[ch] = "C" + ch + "-" + rawTitle;
        adjustChannelStack(channelTitles[ch], stackName, "channel " + ch, useAutoBC);
    }

    Dialog.create("Select channels to keep for " + stackName);
    Dialog.addMessage("Choose which channels to include in the merged composite.\n" +
                      "Channel numbers match C1, C2, ... after splitting.");
    for (ch = 1; ch <= nChannels; ch++)
        Dialog.addCheckbox("Keep channel " + ch, true);
    Dialog.show();

    nKeep = 0;
    channelMask = "";
    for (ch = 1; ch <= nChannels; ch++) {
        keepChannel[ch] = Dialog.getCheckbox();
        if (keepChannel[ch]) {
            nKeep++;
            channelMask = channelMask + "1";
        } else
            channelMask = channelMask + "0";
    }

    if (nKeep < 1)
        exit("At least one channel must be kept for " + stackName + ".");

    for (ch = 1; ch <= nChannels; ch++) {
        if (!keepChannel[ch])
            closeIfOpen(channelTitles[ch]);
    }

    if (nKeep == 1) {
        for (ch = 1; ch <= nChannels; ch++) {
            if (keepChannel[ch]) {
                selectWindow(channelTitles[ch]);
                run("Rename...", "title=[" + stackBase + "]");
            }
        }
        return channelMask;
    }

    mergeCmd = "";
    for (ch = 1; ch <= nChannels; ch++) {
        if (keepChannel[ch])
            mergeCmd = mergeCmd + " c" + ch + "=[" + channelTitles[ch] + "]";
    }
    mergeCmd = mergeCmd + " create keep ignore";
    run("Merge Channels...", mergeCmd);
    run("Rename...", "title=[" + stackBase + "]");

    for (ch = 1; ch <= nChannels; ch++) {
        if (keepChannel[ch])
            closeIfOpen(channelTitles[ch]);
    }

    return channelMask;
}

function adjustChannelStack(winTitle, stackName, channelLabel, useAutoBC) {
    if (!isOpen(winTitle))
        exit("Expected channel window not found: " + winTitle);

    selectWindow(winTitle);
    getDimensions(w, h, c, z, t);
    if (z < 1) z = nSlices;
    if (t < 1) t = 1;

    currentSlice = getSliceNumber();

    if (useAutoBC)
        run("Enhance Contrast", "saturated=0.35 process_all");

    run("Brightness/Contrast...");
    autoNote = "";
    if (useAutoBC)
        autoNote = "Auto brightness/contrast (Enhance Contrast) was applied to all z-planes.\n";
    waitForUser("Review " + channelLabel + " for " + stackName,
        autoNote +
        "Review or fine-tune " + channelLabel + " in the Brightness/Contrast window.\n" +
        "When satisfied, click OK.\n" +
        "The chosen range will be propagated to all z-planes of this channel stack.");

    getMinAndMax(savedMin, savedMax);

    for (zz = 1; zz <= z; zz++) {
        setSlice(zz);
        setMinAndMax(savedMin, savedMax);
    }

    setSlice(currentSlice);
}

function addPersistentMarkerAllZSlices(xpx, ypx, zCount, labelText) {
    getDimensions(w, h, c, z, t);
    if (zCount < 1) zCount = 1;
    if (t < 1) t = 1;

    r = 8;
    if (r < 4) r = 4;

    labelX = xpx + r + 4;
    labelY = ypx - r - 2;
    if (labelY < 12) labelY = ypx + r + 12;

    if (c > 1 || t > 1) {
        Stack.getPosition(oldC, oldZ, oldT);
    } else {
        oldZ = getSliceNumber();
    }

    setFont("SansSerif", 12, "antialiased");

    for (tt = 1; tt <= t; tt++) {
        for (zz = 1; zz <= zCount; zz++) {

            if (c > 1 || t > 1)
                Stack.setPosition(1, zz, tt);
            else
                setSlice(zz);

            makeOval(xpx - r, ypx - r, 2 * r, 2 * r);
            run("Properties...", "stroke=yellow width=2");
            run("Add Selection...");

            makeLine(xpx - r, ypx, xpx + r, ypx);
            run("Properties...", "stroke=yellow width=2");
            run("Add Selection...");

            makeLine(xpx, ypx - r, xpx, ypx + r);
            run("Properties...", "stroke=yellow width=2");
            run("Add Selection...");

            makeText(labelText, labelX, labelY);
            run("Properties...", "stroke=yellow width=1");
            run("Add Selection...");
        }
    }

    run("Select None");
    run("Show Overlay");

    if (c > 1 || t > 1)
        Stack.setPosition(oldC, oldZ, oldT);
    else
        setSlice(oldZ);
}

function closeIfOpen(winTitle) {
    if (isOpen(winTitle)) {
        selectWindow(winTitle);
        close();
    }
}

function getCurrentZSlice() {
    getDimensions(w, h, c, z, t);
    if (c > 1 || t > 1) {
        Stack.getPosition(currC, currZ, currT);
        return currZ;
    } else {
        return getSliceNumber();
    }
}

function setCurrentZSlice(zWanted) {
    getDimensions(w, h, c, z, t);
    if (c > 1 || t > 1) {
        Stack.getPosition(currC, currZ, currT);
        Stack.setPosition(currC, zWanted, currT);
    } else {
        setSlice(zWanted);
    }
}

function sanitizeName(txt) {
    t = txt;
    t = replace(t, " ", "_");
    t = replace(t, "/", "_");
    t = replace(t, "\\", "_");
    t = replace(t, ":", "_");
    t = replace(t, ";", "_");
    t = replace(t, ",", "_");
    t = replace(t, ".", "_");
    t = replace(t, "(", "_");
    t = replace(t, ")", "_");
    t = replace(t, "[", "_");
    t = replace(t, "]", "_");
    t = replace(t, "{", "_");
    t = replace(t, "}", "_");
    t = replace(t, "\"", "_");
    t = replace(t, "'", "_");
    return t;
}

function csvSafe(txt) {
    s = "" + txt;
    s = replace(s, ";", ",");
    s = replace(s, "\n", " ");
    s = replace(s, "\r", " ");
    s = replace(s, "\"", "'");
    return s;
}

function stripExtension(name) {
    dot = lastIndexOf(name, ".");
    if (dot == -1) return name;
    return substring(name, 0, dot);
}

function getSeriesBase(baseName) {
    p = lastIndexOf(baseName, "_z");
    if (p == -1) return baseName;

    token = substring(baseName, p + 2);
    if (lengthOf(token) == 2 && isDigits(token))
        return substring(baseName, 0, p);
    else
        return baseName;
}

function intToStr(v) {
    return "" + floor(v);
}

function numToStr(v) {
    if (v < 0) return "";
    s = d2s(v, numDigits);
    s = replace(s, ".", ",");
    return s;
}

function isLikelyCompanionSlice(fname) {
    base = stripExtension(fname);
    p = lastIndexOf(base, "_z");
    if (p == -1) return 0;

    token = substring(base, p + 2);
    if (lengthOf(token) != 2) return 0;
    if (isDigits(token) == 0) return 0;

    if (token == "00")
        return 0;
    else
        return 1;
}

function isDigits(s) {
    for (kk = 0; kk < lengthOf(s); kk++) {
        ch = substring(s, kk, kk + 1);
        if (indexOf("0123456789", ch) == -1)
            return 0;
    }
    return 1;
}

function absVal(v) {
    if (v < 0) return -v;
    return v;
}

function normalizeUnit(u) {
    s = "" + u;
    s = replace(s, "µ", "u");
    if (s == "um") return "um";
    if (s == "micron") return "um";
    if (s == "microns") return "um";
    if (s == "micrometer") return "um";
    if (s == "micrometers") return "um";
    return s;
}

function parseNumberSmart(s) {
    t = "" + s;
    t = replace(t, " ", "");
    if (indexOf(t, ",") != -1 && indexOf(t, ".") != -1) {
        t = replace(t, ",", "");
    } else if (indexOf(t, ",") != -1) {
        t = replace(t, ",", ".");
    }
    return parseFloat(t);
}

function getDimensionTag(xml, dimID) {
    needle = "DimID=\"" + dimID + "\"";
    p = indexOf(xml, needle);
    if (p == -1) return "";

    prefix = substring(xml, 0, p);
    start = lastIndexOf(prefix, "<DimensionDescription");
    if (start == -1) return "";

    suffix = substring(xml, p);
    endRel = indexOf(suffix, ">");
    if (endRel == -1) return "";

    endPos = p + endRel + 1;
    return substring(xml, start, endPos);
}

function getXmlAttr(tag, attr) {
    patt = attr + "=\"";
    p = indexOf(tag, patt);
    if (p == -1) return "";

    tail = substring(tag, p + lengthOf(patt));
    q = indexOf(tail, "\"");
    if (q == -1) return "";

    return substring(tail, 0, q);
}
