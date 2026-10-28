macro "Measure distances from TIFF and save CSV" {

    requires("1.52p");

    // ----------------------------
    // 1) Select parent folder
    // ----------------------------
    parentDir = getDirectory("Choose the folder containing the TIFF file(s) and metadata folder");
    if (parentDir == "")
        exit("No folder selected.");

    // ----------------------------
    // 2) Find metadata XML
    // ----------------------------
    xmlPath = findPropertiesXML(parentDir);
    if (xmlPath == "")
        exit("Could not find a *_Properties.xml file in the selected folder or its direct subfolders.");

    xmlText = File.openAsString(xmlPath);
    if (xmlText == "")
        exit("Could not read metadata XML:\n" + xmlPath);

    // Read X dimension voxel size and unit
    xBlock = getDimensionBlock(xmlText, "X");
    if (xBlock == "")
        exit("Could not find DimID=\"X\" in metadata XML.");

    pixelSize = getAttributeValue(xBlock, "Voxel");
    unit = getAttributeValue(xBlock, "Unit");

    if (pixelSize == "" || unit == "")
        exit("Could not extract Voxel and/or Unit from X dimension in metadata XML.");

    // Convert comma decimal to dot decimal if needed
    pixelSize = replace(pixelSize, ",", ".");
    px = parseFloat(pixelSize);
    if (isNaN(px) || px <= 0)
        exit("Invalid pixel size read from XML: " + pixelSize);

    // ----------------------------
    // 3) Count TIFF files
    // ----------------------------
    tifCount = countTiffFiles(parentDir);
    if (tifCount == 0)
        exit("No .tif or .tiff files found in the selected folder.");

    // ----------------------------
    // 4) Open image
    //    - one TIFF  -> open directly
    //    - many TIFF -> open as image sequence
    // ----------------------------
    if (tifCount == 1) {
        tifPath = getFirstTiffFile(parentDir);
        if (tifPath == "")
            exit("Could not identify the TIFF file to open.");
        open(tifPath);
    } else {
        // Assumes all TIFFs in the selected folder belong to one stack
        run("Image Sequence...", "open=[" + parentDir + "] sort");
    }

    imgTitle = getTitle();

    // ----------------------------
    // 5) Apply calibration from XML
    // ----------------------------
    // distance=1 pixel corresponds to known=pixelSize in given unit
    run("Set Scale...", "distance=1 known=" + d2s(px, 9) + " unit=[" + unit + "]");

    // ----------------------------
    // 6) Prepare ROI Manager / overlay
    // ----------------------------
    roiManager("reset");
    roiManager("Show All with labels");

    // ----------------------------
    // 7) Repeated measurements
    // ----------------------------
    csv = "Measurement_Number,Distance,Unit,PixelSize,Follicle_Type\n";
    measurementNumber = 1;

    follicleOptions = newArray("primary follicle", "primordial follicle", "not distinguishable");

    showMessage(
        "Instructions",
        "1. Use the straight line tool.\n" +
        "2. Navigate freely between slices/channels if needed.\n" +
        "3. Draw a line.\n" +
        "4. Click OK in the next prompt to record the measurement.\n" +
        "5. After each measurement, select the follicle type.\n\n" +
        "The line length will be computed using the pixel size from the XML metadata."
    );

    keepMeasuring = true;

    while (keepMeasuring) {

        waitForUser("Draw a line with the straight line tool, then click OK to record the measurement.");

        if (selectionType() == -1)
            exit("No selection found. Macro stopped.");

        getSelectionCoordinates(xp, yp);

        if (lengthOf(xp) != 2 || lengthOf(yp) != 2)
            exit("The current selection is not a straight 2-point line. Please use the standard straight line tool.");

        dx = xp[1] - xp[0];
        dy = yp[1] - yp[0];
        distPixels = sqrt(dx*dx + dy*dy);
        distCalibrated = distPixels * px;

        // Add ROI to ROI Manager and keep visible
        roiManager("Add");
        roiCount = roiManager("count");
        roiManager("Select", roiCount - 1);
        roiManager("Rename", "Measure_" + measurementNumber);

        // Ask for qualitative label linked to this measurement
        Dialog.create("Follicle type");
        Dialog.addChoice(
            "What kind of follicle did you just measure the distance to the cortex from?",
            follicleOptions,
            follicleOptions[0]
        );
        Dialog.show();
        follicleType = Dialog.getChoice();

        // Append to CSV
        csv = csv + measurementNumber + "," +
                    d2s(distCalibrated, 6) + "," +
                    unit + "," +
                    d2s(px, 6) + "," +
                    csvEscape(follicleType) + "\n";

        // Ask whether to continue
        keepMeasuring = getBoolean(
            "Recorded measurement #" + measurementNumber +
            "\nDistance = " + d2s(distCalibrated, 6) + " " + unit +
            "\nFollicle type = " + follicleType +
            "\n\nAdd another measurement?"
        );

        measurementNumber = measurementNumber + 1;
    }

    // ----------------------------
    // 8) Save CSV once at the end
    // ----------------------------
    defaultCSV = imgTitle + "_measurements.csv";

    saveDir = getDirectory("Choose the folder where the CSV should be saved");
    if (saveDir == "")
        exit("Measurements were collected, but no save folder was selected.");

    Dialog.create("Save CSV");
    Dialog.addString("Filename:", defaultCSV, 40);
    Dialog.show();
    saveName = Dialog.getString();

    if (saveName == "")
        exit("Measurements were collected, but no CSV filename was provided.");

    if (!endsWith(toLowerCase(saveName), ".csv"))
        saveName = saveName + ".csv";

    savePath = saveDir + saveName;

    File.saveString(csv, savePath);

    showMessage(
        "Done",
        "Measurements saved to:\n" + savePath + "\n\n" +
        "Pixel size used: " + d2s(px, 6) + " " + unit + "/pixel\n" +
        "Number of measurements: " + (measurementNumber - 1)
    );
}


// ======================================================
// Helper functions
// ======================================================

function countTiffFiles(dir) {
    list = getFileList(dir);
    count = 0;

    for (i = 0; i < lengthOf(list); i++) {
        path = dir + list[i];
        if (!File.isDirectory(path)) {
            nameLower = toLowerCase(list[i]);
            if (endsWith(nameLower, ".tif") || endsWith(nameLower, ".tiff"))
                count = count + 1;
        }
    }
    return count;
}

function getFirstTiffFile(dir) {
    list = getFileList(dir);

    for (i = 0; i < lengthOf(list); i++) {
        path = dir + list[i];
        if (!File.isDirectory(path)) {
            nameLower = toLowerCase(list[i]);
            if (endsWith(nameLower, ".tif") || endsWith(nameLower, ".tiff"))
                return path;
        }
    }
    return "";
}

function findPropertiesXML(parentDir) {
    // First: look in parent folder itself
    list = getFileList(parentDir);
    for (i = 0; i < lengthOf(list); i++) {
        path = parentDir + list[i];
        if (!File.isDirectory(path)) {
            nameLower = toLowerCase(list[i]);
            if (endsWith(nameLower, "_properties.xml"))
                return path;
        }
    }

    // Second: look in direct subfolders (e.g. metadata folder)
    for (i = 0; i < lengthOf(list); i++) {
        subPath = parentDir + list[i];
        if (File.isDirectory(subPath)) {
            subList = getFileList(subPath);
            for (j = 0; j < lengthOf(subList); j++) {
                filePath = subPath + subList[j];
                if (!File.isDirectory(filePath)) {
                    nameLower = toLowerCase(subList[j]);
                    if (endsWith(nameLower, "_properties.xml"))
                        return filePath;
                }
            }
        }
    }

    return "";
}

function getDimensionBlock(xmlText, dimID) {
    tag = "DimID=\"" + dimID + "\"";
    start = indexOf(xmlText, tag);
    if (start < 0) return "";

    end = indexOf(xmlText, ">", start);
    if (end < 0) return "";

    return substring(xmlText, start, end + 1);
}

function getAttributeValue(text, attributeName) {
    key = attributeName + "=\"";
    p1 = indexOf(text, key);
    if (p1 < 0) return "";

    p1 = p1 + lengthOf(key);
    p2 = indexOf(text, "\"", p1);
    if (p2 < 0) return "";

    return substring(text, p1, p2);
}

function csvEscape(text) {
    text = replace(text, "\"", "\"\"");
    return "\"" + text + "\"";
}