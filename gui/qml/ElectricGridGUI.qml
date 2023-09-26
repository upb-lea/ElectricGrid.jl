import QtQuick 2.6
import QtQuick.Controls 2.3
import QtQuick.Layouts 1.0
import QtQuick.Window 2.2
import org.julialang 1.0

Window {
    id: mainWindow
    visible: true
    width: 1200
    height: 800
    title: "Node Demo"

    property int sourceIndex: 1
    property int loadIndex: 1
    property int cableIndex: 1
    property bool drawCable: false
    property var drawCableNode
    property bool showMenu: false
    property bool showNodeMenuSource: false
    property bool showNodeMenuLoad: false
    property bool showCableMenu: false
    property int menuOffsetX: 0
    property int menuOffsetY: 0
    property bool sourceButtonHighlighted: false
    property bool loadButtonHighlighted: false
    property bool nodeSourceButtonHighlighted: false
    property bool nodeLoadButtonHighlighted: false
    property var rightClickedNode: null
    property var rightClickedCable: null
    property var cable_length_factor: 50.0

    JuliaSignals {
        signal nodesChanged()
        onNodesChanged: redraw()
    }

    ListModel {
        id: nodes
    }

    ListModel {
        id: cables
    }

    function addSource() {
        if (sourceButton.text === "Add Source" && sourceButtonHighlighted) {
            var uid = "source" + sourceIndex;
            sourceIndex += 1;

            var ex = (menuOffsetX + menu.width / 2) / mainWindow.width;
            var ey = (menuOffsetY + menu.height / 2) / mainWindow.height;

            // Generate random pwr value
            var pwr = Math.random() * 100; // Modify the range as needed

            nodes.append({
                "uid": uid,
                "ex": ex,
                "ey": ey,
                "color": "red",
                "source_type": "ideal",
                "control_type": "Classic",
                "mode": "Droop",
                "filter": "L",
                "pwr": pwr,
            });

            Julia.addSource(uid, "ideal", "Classic", "Droop", "L", pwr);
            showMenu = false;
            showNodeMenuSource = false;
        }
    }


    // Updated addLoad() function
    function addLoad() {
        if (loadButton.text === "Add Load" && loadButtonHighlighted) {
            var uid = "load" + loadIndex;
            loadIndex += 1;

            var ex = (menuOffsetX + menu.width / 2) / mainWindow.width;
            var ey = (menuOffsetY + menu.height / 2) / mainWindow.height;

            // Get the selected impedance type index from the ComboBox
            //var impedance_index = impedanceComboBoxLoad.currentIndex;
            // Get the selected impedance type text using the index
            //var impedance = impedanceComboBoxLoad.model.get(impedance_index).impedance;
            // Generate random impedance, mode, pwr, and filter values
            var pwr = Math.random() * 10; // Modify the range as needed

            nodes.append({
                "uid": uid,
                "ex": ex,
                "ey": ey,
                "color": "blue",
                "impedance": "R",
                "pwr": pwr,
                //"loadImpedanceIndex": impedanceComboBoxLoad.currentIndex  // Store impedance index for each load node
                "source_type": "", //dummy entry (needed for icons)
                "control_type": "", //dummy entry (needed for icons)
            });

            Julia.addLoad(uid, "R", pwr);
            showMenu = false;
            showNodeMenuLoad = false;
        }
    }

    function addCable(to_uid) {
        var valid = true;

        if (drawCableNode.uid === to_uid) {
            valid = false;
            return;
        }

        for (var i = 0; i < cables.count; i++) {
            if ((cables.get(i)["from"] === drawCableNode.uid && cables.get(i)["to"] === to_uid) ||
                (cables.get(i)["from"] === to_uid && cables.get(i)["to"] === drawCableNode.uid)) {
                valid = false;
                break;
            }
        }

        if (valid) {
            var uid = "cable" + cableIndex;
            cableIndex += 1;

            var fromNode, toNode;

            // Find the fromNode and toNode based on their uids
            for (var j = 0; j < nodes.count; j++) {
                if (nodes.get(j)["uid"] === drawCableNode.uid) {
                    fromNode = nodes.get(j);
                }
                if (nodes.get(j)["uid"] === to_uid) {
                    toNode = nodes.get(j);
                }
            }

            // Calculate the cable coordinates
            var fromX = fromNode.ex * mainWindow.width + 10;
            var fromY = fromNode.ey * mainWindow.height + 10;
            var toX = toNode.ex * mainWindow.width + 10;
            var toY = toNode.ey * mainWindow.height + 10;

            var length = Math.hypot(fromX - toX, fromY - toY )

            var capacity = parseFloat(param1FieldCable.text);
            var inductance = parseFloat(param2FieldCable.text);
            var resistance = parseFloat(param3FieldCable.text);

            cables.append({
                "uid": uid,
                "from": drawCableNode.uid,
                "to": to_uid,
                "fromX": fromX,
                "fromY": fromY,
                "toX": toX,
                "toY": toY,
                "color": "black",
                "lineWidth": 2,
                "length": length / cable_length_factor,
                "showlength": false,
                "capacity": capacity,
                "inductance": inductance,
                "resistance": resistance
            });

            Julia.addCable(uid, drawCableNode.uid, to_uid, length / cable_length_factor, capacity, inductance, resistance);
            canvas.requestPaint();
        } else {
            // Cancel the ongoing connection
            drawCable = false;
            canvas.requestPaint();
        }
    }

    function updateSource() {
        Julia.updateSource(rightClickedNode.uid,
                            rightClickedNode.source_type,
                            rightClickedNode.control_type,
                            rightClickedNode.mode,
                            rightClickedNode.filter,
                            rightClickedNode.pwr);
    }

    function updateLoad() {
        Julia.updateLoad(rightClickedNode.uid,
                            rightClickedNode.impedance,
                            rightClickedNode.pwr);
    }

    function updateCable() {
        Julia.updateCable(rightClickedCable.uid,
                            rightClickedCable.from,
                            rightClickedCable.to,
                            rightClickedCable.length,
                            rightClickedCable.capacity,
                            rightClickedCable.inductance,
                            rightClickedCable.resistance);
    }

    function deleteRightClickedNode() {
        showNodeMenuLoad = false;
        showNodeMenuSource = false;
        showCableMenu = false;

        var nodeUID = rightClickedNode.uid
        var foundIndex = -1;
        for (var i = 0; i < nodes.count; ++i) {
            if (nodes.get(i).uid === nodeUID) {
                foundIndex = i;
                break;
            }
        }
        nodes.remove(foundIndex)

        var foundIndices = []
        for (var i = 0; i < cables.count; ++i) {
            if (cables.get(i).from === nodeUID || cables.get(i).to === nodeUID) {
                foundIndices.push(i);
            }
        }

        for (var i = foundIndices.length - 1; i >= 0 ; i--) {
            cables.remove(foundIndices[i])
        }

        canvas.requestPaint();

        Julia.deleteNode(nodeUID)
    }

    function deleteRightClickedCable() {
        showNodeMenuLoad = false;
        showNodeMenuSource = false;
        showCableMenu = false;

        var cableUID = rightClickedCable.uid
        var foundIndex = -1;
        for (var i = 0; i < cables.count; ++i) {
            if (cables.get(i).uid === cableUID) {
                foundIndex = i;
                break;
            }
        }
        cables.remove(foundIndex)

        canvas.requestPaint();

        Julia.deleteCable(cableUID)
    }

    Rectangle {
        id: mouseFollower
        width: 0
        height: 0
    }

    Item {
        id: root
        anchors.fill: parent

        MouseArea {
            id: globalMouseArea
            anchors.fill: parent
            hoverEnabled: true
            acceptedButtons: Qt.LeftButton | Qt.RightButton

            onClicked: {
                var clickedNode = isNodeClicked(mouse);

                if (showMenu && !menu.containsMouse) {
                    // Clicked outside the menu, hide the menu
                    showMenu = false;
                } else if (!drawCable && mouse.button === Qt.RightButton && !clickedNode) {
                    // Right-clicked on the canvas, show the menu
                    showMenu = true;
                    menuOffsetX = mouse.x - menu.width / 2; // Set the menu's X position offset
                    menuOffsetY = mouse.y - menu.height / 2; // Set the menu's Y position offset
                } else if (!drawCable && mouse.button === Qt.RightButton && clickedNode) {
                    forceActiveFocus()
                    showNodeMenuLoad = false;
                    showNodeMenuSource = false;
                    showCableMenu = false;

                    //We have to set RightClickedNode to null first to prevent errors
                    rightClickedNode = null;
                    // Right-clicked on the node, show the node menu
                    rightClickedNode = clickedNode;

                    if (rightClickedNode.color === "red") {
                        showNodeMenuSource = true;
                    } else if (rightClickedNode.color === "blue") {
                        showNodeMenuLoad = true;
                    }
                    nodeMenuSource.x = clickedNode.ex * mainWindow.width + 25; 
                    nodeMenuSource.y = clickedNode.ey * mainWindow.height + 10; 
                    nodeMenuLoad.x = clickedNode.ex * mainWindow.width + 25; 
                    nodeMenuLoad.y = clickedNode.ey * mainWindow.height + 10; 
                } else if (drawCable) {
                    // Clicked anywhere on the window other than a node, cancel the ongoing connection
                    drawCable = false;
                    canvas.requestPaint();
                }

                if (showNodeMenuSource && !nodeMenuSource.containsMouse && !clickedNode) {
                    // Left-clicked outside the Node menu for source nodes, hide the node menu
                    forceActiveFocus()
                    showNodeMenuSource = false;
                    rightClickedNode = null;
                    canvas.requestPaint();
                }

                if (showNodeMenuLoad && !nodeMenuLoad.containsMouse && !clickedNode) {
                    // Left-clicked outside the Node menu for load nodes, hide the node menu
                    forceActiveFocus()
                    showNodeMenuLoad = false;
                    rightClickedNode = null;
                    canvas.requestPaint();
                }
                
                if (showCableMenu && !showCableMenu.containsMouse && !clickedNode) {
                    // Left-clicked outside the Node menu for load nodes, hide the node menu
                    forceActiveFocus()
                    showCableMenu = false;
                    rightClickedCable = null;
                    canvas.requestPaint();
                }
            }

            function isNodeClicked(mouse) {
                for (var i = 0; i < nodes.count; i++) {
                    var node = nodes.get(i);
                    var nodeRect = getNodeRect(node);
                    if (
                        mouse.x >= nodeRect.x &&
                        mouse.x <= nodeRect.x + nodeRect.width &&
                        mouse.y >= nodeRect.y &&
                        mouse.y <= nodeRect.y + nodeRect.height
                    ) {
                        return nodes.get(i); // Return the node data instead of true
                    }
                }
                return null; // Return null if no node is clicked
            }

            function getNodeRect(node) {
                var nodeX = node.ex * mainWindow.width;
                var nodeY = node.ey * mainWindow.height;
                var nodeSize = 22;
                return Qt.rect(nodeX, nodeY, nodeSize, nodeSize);
            }

            onPositionChanged: {
                if (drawCable) {
                    mouseFollower.x = mouse.x;
                    mouseFollower.y = mouse.y;
                    canvas.requestPaint();
                }
            }
        }
    }

    function getConnectedNodes(uid) {
        var connectedNodes = [];
        for (var i = 0; i < cables.count; i++) {
            if (cables.get(i).from === uid) {
                connectedNodes.push(cables.get(i).to);
            } else if (cables.get(i).to === uid) {
                connectedNodes.push(cables.get(i).from);
            }
        }
        return connectedNodes;
    }


    Repeater {
        anchors.fill: parent
        model: cables

        Rectangle {
            x: model.fromX
            y: model.fromY - 3

            property var len: Math.hypot(model.fromX - model.toX, model.fromY - model.toY )

            width: len
            height: 6
            color: "transparent"
            transform: Rotation {origin.x: 0; origin.y: 3; angle: Math.atan2(model.toY - model.fromY, model.toX - model.fromX) * 180 / Math.PI}

            MouseArea {
                anchors.fill: parent
                hoverEnabled: true

                onEntered: parent.color = "black"
                onExited: parent.color = "transparent"

                acceptedButtons: Qt.RightButton | Qt.LeftButton

                onClicked: {
                    showNodeMenuLoad = false;
                    showNodeMenuSource = false;
                    if (mouse.button === Qt.RightButton) {
                        // Show the cable menu if right-clicked on the cable
                        showCableMenu = true;
                        rightClickedCable = model; // Store the clicked cable
                    } else if (showCableMenu && !cableMenuArea.containsMouse) {
                        // Left-clicked outside the cable menu, hide the cable menu
                        showCableMenu = false;
                        rightClickedCable = null;
                        canvas.requestPaint();
                    }
                }
            }

            Text {
                x: model.length * cable_length_factor / 2
                y: 5
                rotation: (model.toX - model.fromX) < 0 ? 180 : 0
                text: model.length.toFixed(2)
                visible: model.showlength
                font.bold: true
                verticalAlignment: Text.AlignVCenter
                Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
            }
        }
    }

    Canvas {
        id: canvas
        anchors.fill: parent

        onPaint: {
            var ctx = getContext("2d");
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            for (var i = 0; i < cables.count; i++) {
                ctx.beginPath();
                ctx.strokeStyle = cables.get(i)["color"]; // Use the color from the cables array
                ctx.lineWidth = cables.get(i)["lineWidth"];
                
                var fromX = cables.get(i)["fromX"];
                var fromY = cables.get(i)["fromY"];
                var toX = cables.get(i)["toX"];
                var toY = cables.get(i)["toY"];

                ctx.moveTo(fromX, fromY);
                ctx.lineTo(toX, toY);
                ctx.stroke();
            }

            if (drawCable) {
                // Draw the black cable on top of the transparent cable
                ctx.beginPath();
                ctx.strokeStyle = "black";
                ctx.lineWidth = 2;
                
                ctx.moveTo(drawCableNode["x"], drawCableNode["y"]);
                ctx.lineTo(mouseFollower.x, mouseFollower.y);
                ctx.stroke();
            }
        }
    }


    function redraw() {
        canvas.requestPaint();
    }


    Repeater {
        anchors.fill: parent
        model: nodes

        Rectangle {
            width: 20
            height: 20
            color: model.color
            radius: width / 2

            x: model.ex * mainWindow.width
            y: model.ey * mainWindow.height

            Image {
                x: model.color === "red" ? -15 : -5
                y: -40
                width: model.color === "red" ? 25 : 30
                height: model.color === "red" ? 25 : 30
                source: model.color === "red" ? (model.source_type === "ideal" ? "ideal.png" : (model.source_type === "pv" ? "pv.png" : "battery.png")) : "load.png"

                MouseArea { 
                    anchors.fill:parent;
                    onClicked: {
                        print(model.color === "red" ? (model.source_type === "ideal" ? "ideal.png" : (model.source_type === "pv" ? "pv.png" : "battery.png")) : "load.png")
                    }
                }
            }

            Image {
                x: 15
                y: -40
                width: 25
                height: 25
                source: model.color === "red" ? (model.control_type === "Classic" ? "classic.png" : "robot-solid.png") : "load.png"
                visible: model.color === "red"

                MouseArea { 
                    anchors.fill:parent;
                    onClicked: {
                        print(model.color === "red" ? (model.control_type) : "load.png")
                    }
                }

            }

            MouseArea {
                anchors.fill: parent
                drag.target: parent

                onClicked: {
                    if (drawCable) {
                        drawCable = false;
                        addCable(model.uid);
                    } else {
                        drawCableNode = {
                            "uid": model.uid,
                            "x": parent.x + 10,
                            "y": parent.y + 10
                        };
                        drawCable = true;
                    }
                    //addLoad(impedanceComboBoxLoad.model.get(impedanceComboBoxLoad.currentIndex).impedance);
                }
                onPositionChanged: {
                    drawCable = false;
                    model.ex = parent.x / mainWindow.width;
                    model.ey = parent.y / mainWindow.height;

                    // Update the cable coordinates
                    for (var i = 0; i < cables.count; i++) {
                        if (cables.get(i)["from"] === model.uid || cables.get(i)["to"] === model.uid) {
                            var fromX, fromY, toX, toY;

                            for (var j = 0; j < nodes.count; j++) {
                                if (nodes.get(j)["uid"] === cables.get(i)["from"]) {
                                    fromX = nodes.get(j)["ex"] * mainWindow.width + 10;
                                    fromY = nodes.get(j)["ey"] * mainWindow.height + 10;
                                }
                                if (nodes.get(j)["uid"] === cables.get(i)["to"]) {
                                    toX = nodes.get(j)["ex"] * mainWindow.width + 10;
                                    toY = nodes.get(j)["ey"] * mainWindow.height + 10;
                                }
                            }

                            var length = Math.hypot(fromX - toX, fromY - toY )

                            cables.setProperty(i, "fromX", fromX);
                            cables.setProperty(i, "fromY", fromY);
                            cables.setProperty(i, "toX", toX);
                            cables.setProperty(i, "toY", toY);
                            cables.setProperty(i, "length", length / cable_length_factor);
                        }
                    }

                    // Update the node menu position
                    if (rightClickedNode && rightClickedNode.uid === model.uid) {
                        nodeMenuSource.x = parent.x + 25; // Adjust the x position offset to align the menu to the right of the node
                        nodeMenuSource.y = parent.y + 10; // Adjust the y position offset to center the menu vertically

                        nodeMenuLoad.x = parent.x + 25; // Adjust the x position offset to align the menu to the right of the node
                        nodeMenuLoad.y = parent.y + 10; // Adjust the y position offset to center the menu vertically
                    }

                    redraw(); // Request a redraw to update the canvas
                }
                onPressed:{
                    for (var i = 0; i < cables.count; i++) {
                        if (cables.get(i)["from"] === model.uid || cables.get(i)["to"] === model.uid) {
                            cables.setProperty(i, "showlength", true);
                        }
                    }
                }
                onReleased:{
                    for (var i = 0; i < cables.count; i++) {
                        if (cables.get(i)["from"] === model.uid || cables.get(i)["to"] === model.uid) {
                            cables.setProperty(i, "showlength", false);

                            Julia.updateCable(cables.get(i)["uid"],
                                                cables.get(i)["from"],
                                                cables.get(i)["to"],
                                                cables.get(i)["length"],
                                                cables.get(i)["capacity"],
                                                cables.get(i)["inductance"],
                                                cables.get(i)["resistance"]);
                        }
                    }
                }
            }
        }
    }


    Rectangle {
        id: menu
        width: 100
        height: 100
        color: "lightgray"
        radius: 10
        visible: showMenu

        // Set the position of the menu based on the menuOffsetX and menuOffsetY properties
        x: mainWindow.menuOffsetX
        y: mainWindow.menuOffsetY

        ColumnLayout {
            anchors.fill: parent
            spacing: 0 // Remove the spacing between the buttons

            Button {
                id: sourceButton
                Layout.fillWidth: true
                text: "Add Source"
                onClicked: mainWindow.addSource()

                // Add MouseArea for hover detection
                MouseArea {
                    anchors.fill: parent
                    hoverEnabled: true

                    // Change the background color when hovered over
                    onEntered: {
                        sourceButton.background.color = "lightblue"
                        mainWindow.sourceButtonHighlighted = true
                    }
                    onExited: {
                        sourceButton.background.color = "transparent"
                        mainWindow.sourceButtonHighlighted = false
                    }
                }

                background: Rectangle {
                    color: "transparent"
                }
            }

            Button {
                id: loadButton
                Layout.fillWidth: true
                text: "Add Load"
                onClicked: mainWindow.addLoad()

                // Add MouseArea for hover detection
                MouseArea {
                    anchors.fill: parent
                    hoverEnabled: true

                    // Change the background color when hovered over
                    onEntered: {
                        loadButton.background.color = "lightblue"
                        mainWindow.loadButtonHighlighted = true
                    }
                    onExited: {
                        loadButton.background.color = "transparent"
                        mainWindow.loadButtonHighlighted = false
                    }
                }

                background: Rectangle {
                    color: "transparent"
                }
            }
        }

        // MouseArea for "Add Source" button
        MouseArea {
            width: menu.width
            height: menu.height / 2
            anchors.top: menu.top
            onClicked: mainWindow.addSource()
        }

        // MouseArea for "Add Load" button
        MouseArea {
            width: menu.width
            height: menu.height / 2
            anchors.bottom: menu.bottom
            onClicked: mainWindow.addLoad()
        }
    }

    Rectangle {
        id: nodeMenuSource
        width: 250
        height: (sourceTypeComboBoxSource.height + controlTypeComboBoxSource.height + (modeText.visible && modeComboBoxSource.visible ? modeComboBoxSource.height : 0) + param3FieldSource.height + filterComboBoxSource.height + 100) // Adjust the height based on the text field heights        
        color: "lightgray"
        radius: 10
        visible: showNodeMenuSource

        // Set the position of the menu relative to the clicked node
        x: rightClickedNode ? rightClickedNode.x + 10 : 0
        y: rightClickedNode ? rightClickedNode.y + 10 : 0

        property var control_type: "Classic" // Set default value to "Classic"
        property var mode: "Droop" // Set default value to "Droop"
        property var pwr: 0
        property var filter: "L" // Set default value to "L"

        

        ColumnLayout {
            id: contentColumnSource
            anchors.fill: parent
            spacing: 10 // Add spacing between the items

            Item {
                height: 5 // Adjust the top spacing
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Source_Type:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                // Dropdown ComboBox for the Source Type
                ComboBox {
                    id: sourceTypeComboBoxSource
                    currentIndex: rightClickedNode ? (rightClickedNode.source_type ? 
                    (rightClickedNode.source_type === "ideal" ? 0 : 
                    (rightClickedNode.source_type === "pv" ? 1 : 2)) : 0) : 0

                    Layout.fillWidth: true
                    model: ListModel {
                        ListElement {
                            sourceType: "ideal"
                        }
                        ListElement {
                            sourceType: "pv"
                        }
                        ListElement {
                            sourceType: "battery"
                        }
                    }
                    textRole: "sourceType" // Set the textRole to the property name

                    onCurrentIndexChanged: {
                        if(rightClickedNode) {
                            rightClickedNode.source_type = sourceTypeComboBoxSource.model.get(sourceTypeComboBoxSource.currentIndex).sourceType;

                            updateSource();
                        } 
                    }
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Ctrl_Type:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                // Dropdown ComboBox for the Control Type
                ComboBox {
                    id: controlTypeComboBoxSource
                    currentIndex: rightClickedNode ? (rightClickedNode.control_type ? 
                    (rightClickedNode.control_type === "Classic" ? 0 : 1) : 0) : 0

                    Layout.fillWidth: true
                    model: ListModel {
                        ListElement {
                            controlType: "Classic"
                        }
                        ListElement {
                            controlType: "RL"
                        }
                    }
                    textRole: "controlType" // Set the textRole to the property name

                    onCurrentIndexChanged: {
                        if(rightClickedNode) {
                            rightClickedNode.control_type = controlTypeComboBoxSource.model.get(controlTypeComboBoxSource.currentIndex).controlType;

                            // Hide the mode text and ComboBox when control_type is "RL"
                            if (rightClickedNode.control_type === "RL") {
                                modeText.visible = false;
                                modeComboBoxSource.visible = false;
                            } else {
                                modeText.visible = true;
                                modeComboBoxSource.visible = true;
                            }

                            updateSource();
                        } else {
                            modeText.visible = true;
                            modeComboBoxSource.visible = true;
                        }
                    }
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    id: modeText
                    text: "Mode:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                // Dropdown ComboBox for the Mode
                ComboBox {
                    id: modeComboBoxSource
                    currentIndex: rightClickedNode ? (rightClickedNode.mode ? 
                    (rightClickedNode.mode === "Droop" ? 0 : 
                    (rightClickedNode.mode === "PQ" ? 1 : 
                    (rightClickedNode.mode === "Swing" ? 2 : 3))) : 0) : 0

                    Layout.fillWidth: true
                    model: ListModel {
                        ListElement {
                            mode: "Droop"
                        }
                        ListElement {
                            mode: "PQ"
                        }
                        ListElement {
                            mode: "Swing"
                        }
                        ListElement {
                            mode: "Synchronverter"
                        }
                    }
                    textRole: "mode" // Set the textRole to the property name

                    onCurrentIndexChanged: {
                        if(rightClickedNode) {
                            rightClickedNode.mode = modeComboBoxSource.model.get(modeComboBoxSource.currentIndex).mode;

                            updateSource();
                        }
                    }
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Pwr:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                MouseArea {
                    Layout.fillWidth: true // Fill the available space
                    Layout.preferredHeight: param3FieldSource.height // Set the preferred height
                    preventStealing: true // Prevent the MouseArea from stealing focus
                    onClicked: {
                        if (!param3FieldSource.activeFocus) {
                            param3FieldSource.forceActiveFocus();
                        }
                    }
                }

                TextField {
                    id: param3FieldSource
                    text: rightClickedNode ? (rightClickedNode.pwr ? rightClickedNode.pwr.toFixed(2) : "0") :  "0"
                    Layout.fillWidth: true
                    placeholderText: "random"
                    focus: true // Set focus on the text field
                    readOnly: false // Make the text field editable
                    onAccepted: {
                        // Handle the accepted event
                        // Update the stored value when editing is finished
                        if (rightClickedNode) {
                            rightClickedNode.pwr = parseFloat(text);

                            updateSource();
                        }
                    }
                    onEditingFinished: {
                        // Handle the Enter key pressed event
                        // Update the stored value when editing is finished
                        if (rightClickedNode) {
                            rightClickedNode.pwr = parseFloat(text);

                            updateSource();
                        }
                    }
                }

                Text {
                    text: " VA"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Filter:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                // Dropdown ComboBox for the Filter
                ComboBox {
                    id: filterComboBoxSource

                    currentIndex: rightClickedNode ? (rightClickedNode.filter ? 
                    (rightClickedNode.filter === "L" ? 0 : 
                    (rightClickedNode.filter === "LC" ? 1 : 2)) : 0) : 0

                    Layout.fillWidth: true
                    model: ListModel {
                        ListElement {
                            filter: "L"
                        }
                        ListElement {
                            filter: "LC"
                        }
                        ListElement {
                            filter: "LCL"
                        }
                    }
                    textRole: "filter" // Set the textRole to the property name

                    onCurrentIndexChanged: {
                        if(rightClickedNode) {
                            rightClickedNode.filter = filterComboBoxSource.model.get(filterComboBoxSource.currentIndex).filter;

                            updateSource();
                        }
                    }

                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                        width: 5 // Spacer width
                    }

                    Image {
                    width: 20
                    height: 20
                    source: "trash-can-solid-small.png"
                    MouseArea { anchors.fill:parent;
                        onClicked: {
                            deleteRightClickedNode()
                        }
                    }
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            Item {
                height: 15
            }
        }

        // Adjust the right edge spacing
        anchors.rightMargin: 5
    }

    Rectangle {
        id: nodeMenuLoad
        width: 200
        height: impedanceComboBoxLoad.height + param3FieldLoad.height + 80 // Adjust the height based on the text field heights
        color: "lightgray"
        radius: 10
        visible: showNodeMenuLoad

        // Set the position of the menu relative to the clicked node
        x: rightClickedNode ? rightClickedNode.x + 10 : 0
        y: rightClickedNode ? rightClickedNode.y + 10 : 0

        ColumnLayout {
            id: contentColumnLoad
            anchors.fill: parent
            spacing: 10 // Add spacing between the items

            Item {
                height: 5 // Adjust the top spacing
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Impedance:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                // Dropdown ComboBox for the Control Type
                ComboBox {
                    id: impedanceComboBoxLoad

                    currentIndex: rightClickedNode ? (rightClickedNode.impedance ? 
                    (rightClickedNode.impedance === "R" ? 0 : 
                    (rightClickedNode.impedance === "L" ? 1 : 
                    (rightClickedNode.impedance === "C" ? 2 : 
                    (rightClickedNode.impedance === "LC" ? 3 : 
                    (rightClickedNode.impedance === "RC" ? 4 : 
                    (rightClickedNode.impedance === "RL" ? 5 : 6)))))) : 6) : 6

                    Layout.fillWidth: true
                    model: ListModel {
                        ListElement { impedance: "R" }
                        ListElement { impedance: "L" }
                        ListElement { impedance: "C" }
                        ListElement { impedance: "LC" }
                        ListElement { impedance: "RC" }
                        ListElement { impedance: "RL" }
                        ListElement { impedance: "RLC" }
                    }
                    textRole: "impedance" // Set the textRole to the property name

                    onCurrentIndexChanged: {
                        if(rightClickedNode) {
                            rightClickedNode.impedance = impedanceComboBoxLoad.model.get(impedanceComboBoxLoad.currentIndex).impedance;

                            updateLoad();
                        }
                    }
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Pwr:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                MouseArea {
                    Layout.fillWidth: true // Fill the available space
                    Layout.preferredHeight: param3FieldLoad.height // Set the preferred height
                    preventStealing: true // Prevent the MouseArea from stealing focus
                    onClicked: {
                        if (!param3FieldLoad.activeFocus) {
                            param3FieldLoad.forceActiveFocus();
                        }
                    }
                }

                TextField {
                    id: param3FieldLoad
                    text: rightClickedNode ? (rightClickedNode.pwr ? rightClickedNode.pwr.toFixed(2) : "0") :  "0" // Use stored value as default or set to "0" if undefined
                    Layout.fillWidth: true
                    placeholderText: "Pwr"
                    focus: true // Set focus on the text field
                    readOnly: false // Make the text field editable
                    onAccepted: {
                        // Handle the accepted event
                        // Update the stored value when editing is finished
                        if (rightClickedNode) {
                            rightClickedNode.pwr = parseFloat(text);

                            updateLoad();
                        }
                    }
                    onEditingFinished: {
                        // Handle the Enter key pressed event
                        // Update the stored value when editing is finished
                        if (rightClickedNode) {
                            rightClickedNode.pwr = parseFloat(text);

                            updateLoad();
                        }
                    }
                }

                Text {
                    text: " VA"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                        width: 5 // Spacer width
                    }

                    Image {
                    width: 20
                    height: 20
                    source: "trash-can-solid-small.png"
                    MouseArea { anchors.fill:parent;
                        onClicked: {
                            deleteRightClickedNode()
                        }
                    }
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            Item {
                height: 15
            }
        }
        // Adjust the right edge spacing
        anchors.rightMargin: 5
    }   

    Rectangle {
        id: cableMenu
        width: 300
        height: lengthFieldCable.height + param1FieldCable.height + param2FieldCable.height + param3FieldCable.height + 100 // Adjust the height based on the text field heights
        color: "lightgray"
        radius: 10
        visible: showCableMenu

        // Set the position of the menu relative to the clicked cable
        x: rightClickedCable ? (rightClickedCable.fromX + rightClickedCable.toX) / 2 + 10 : 0
        y: rightClickedCable ? (rightClickedCable.fromY + rightClickedCable.toY) / 2 + 10 : 0

        MouseArea {
            anchors.fill: parent
            onClicked: {
                forceActiveFocus()
            }
        }

        ColumnLayout {
            id: contentColumnCable
            anchors.fill: parent
            spacing: 10 // Add spacing between the items

            Item {
                height: 5 // Adjust the top spacing
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Length:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                Text {
                    id: lengthFieldCable
                    text: (rightClickedCable ? rightClickedCable.length.toFixed(2) :  "1.0") + " km"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Capacity Coatings:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                MouseArea {
                    Layout.fillWidth: true // Fill the available space
                    Layout.preferredHeight: param1FieldCable.height // Set the preferred height
                    preventStealing: true // Prevent the MouseArea from stealing focus
                    onClicked: {
                        if (!param1FieldCable.activeFocus) {
                            param1FieldCable.forceActiveFocus();
                        }
                    }
                }

                TextField {
                    id: param1FieldCable
                    text: rightClickedCable ? (rightClickedCable.capacity ? rightClickedCable.capacity.toFixed(2) : "0.4") :  "0.4" // Use stored value as default or set to "0" if undefined
                    Layout.fillWidth: true
                    placeholderText: "Capacity"
                    focus: true // Set focus on the text field
                    readOnly: false // Make the text field editable
                    onAccepted: {
                        // Handle the accepted event
                        // Update the stored value when editing is finished
                        if (rightClickedCable) {
                            rightClickedCable.capacity = parseFloat(text);

                            updateCable()
                        }
                    }
                    onEditingFinished: {
                        // Handle the Enter key pressed event
                        // Update the stored value when editing is finished
                        if (rightClickedCable) {
                            rightClickedCable.capacity = parseFloat(text);

                            updateCable()
                        }
                    }
                }

                Text {
                    id: capacityFieldCable
                    text: " -> C: " + (rightClickedCable ? (rightClickedCable.length * rightClickedCable.capacity).toFixed(2) :  "1.0") + " µF"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "Operating Inductor:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                MouseArea {
                    Layout.fillWidth: true // Fill the available space
                    Layout.preferredHeight: param2FieldCable.height // Set the preferred height
                    preventStealing: true // Prevent the MouseArea from stealing focus
                    onClicked: {
                        if (!param2FieldCable.activeFocus) {
                            param2FieldCable.forceActiveFocus();
                        }
                    }
                }

                TextField {
                    id: param2FieldCable
                    text: rightClickedCable ? (rightClickedCable.inductance ? rightClickedCable.inductance.toFixed(2) : "0.264") :  "0.264" // Use stored value as default or set to "0" if undefined
                    Layout.fillWidth: true
                    placeholderText: "Inductance"
                    focus: true // Set focus on the text field
                    readOnly: false // Make the text field editable
                    onAccepted: {
                        // Handle the accepted event
                        // Update the stored value when editing is finished
                        if (rightClickedCable) {
                            rightClickedCable.inductance = parseFloat(text);

                            updateCable()
                        }
                    }
                    onEditingFinished: {
                        // Handle the Enter key pressed event
                        // Update the stored value when editing is finished
                        if (rightClickedCable) {
                            rightClickedCable.inductance = parseFloat(text);

                            updateCable()
                        }
                    }
                }

                Text {
                    id: inductanceFieldCable
                    text: " -> L: " + (rightClickedCable ? (rightClickedCable.length * rightClickedCable.inductance).toFixed(2) :  "1.0") + " mH"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                    width: 5 // Spacer width
                }

                Text {
                    text: "AC Resistor:"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                MouseArea {
                    Layout.fillWidth: true // Fill the available space
                    Layout.preferredHeight: param3FieldCable.height // Set the preferred height
                    preventStealing: true // Prevent the MouseArea from stealing focus
                    onClicked: {
                        if (!param3FieldCable.activeFocus) {
                            param3FieldCable.forceActiveFocus();
                        }
                    }
                }

                TextField {
                    id: param3FieldCable
                    text: rightClickedCable ? (rightClickedCable.resistance ? rightClickedCable.resistance.toFixed(2) : "0.722") :  "0.722" // Use stored value as default or set to "0" if undefined
                    Layout.fillWidth: true
                    placeholderText: "Resistance"
                    focus: true // Set focus on the text field
                    readOnly: false // Make the text field editable
                    onAccepted: {
                        // Handle the accepted event
                        // Update the stored value when editing is finished
                        if (rightClickedCable) {
                            rightClickedCable.resistance = parseFloat(text);

                            updateCable()
                        }
                    }
                    onEditingFinished: {
                        // Handle the Enter key pressed event
                        // Update the stored value when editing is finished
                        if (rightClickedCable) {
                            rightClickedCable.resistance = parseFloat(text);

                            updateCable()
                        }
                    }
                }

                Text {
                    id: resistanceFieldCable
                    text: " -> R: " + (rightClickedCable ? (rightClickedCable.length * rightClickedCable.resistance).toFixed(2) :  "1.0") + " Ohm"
                    font.bold: true
                    verticalAlignment: Text.AlignVCenter
                    Layout.alignment: Qt.AlignLeft // Use Layout.alignment instead of anchors
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            RowLayout {
                spacing: 10

                Item {
                        width: 5 // Spacer width
                    }

                    Image {
                    width: 20
                    height: 20
                    source: "trash-can-solid-small.png"
                    MouseArea { anchors.fill:parent;
                        onClicked: {
                            deleteRightClickedCable()
                        }
                    }
                }

                Item {
                    width: 5 // Spacer width
                }
            }

            Item {
                height: 15
            }

        }
        // Adjust the right edge spacing
        anchors.rightMargin: 5
    }     
}

