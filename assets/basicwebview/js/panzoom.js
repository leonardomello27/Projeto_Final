// Pan-Zoom for model view

// Copyright 2021 MathWorks, Inc.

// ========================================================================

const modelStage = document.querySelector('.model-stage');
const modelImage = document.querySelector('.model-image');
const imageWrap = document.querySelector('.image-wrap');
const resetButton = document.querySelector('.reset');
const zoomInButton = document.querySelector('.zoom-in');
const zoomOutButton = document.querySelector('.zoom-out');

let offsetX = 0;
let offsetY = 0;

let posX = 0;
let posY = 0;

let startPanX = 0;
let startPanY = 0;

let isDown = false;
let scale = 1;

let zoomInBtnPressed = false;
let zoomOutBtnPressed = false;

let zoomInTimer;
let zoomOutTimer;

let zeroOriginActive = false;
let switchingToWheelScroll = false;

let justLoaded = true;

let lastScrollPosX = 0;
let lastScrollPosY = 0;

// ========================================================================

function setup(htmlComponent) {
    htmlComponent.addEventListener("DataChanged", function(event) {
        resetView();
    });     
}

// ========================================================================

function mousePressed(event) {
    captureMousePosition(event);
    startPanX = posX;
    startPanY = posY;
    isDown = true;
}

// ========================================================================

function disable(event) {
    isDown = false;
}

// ========================================================================

function pan(event) {
    event.preventDefault();

    if(!isDown) return;

    if (event.clientX >= window.innerWidth - 100 &&
        event.clientY >= window.innerHeight - 50) return;

    offsetX += event.clientX - startPanX;
    offsetY += event.clientY - startPanY;

    modelStage.style.transform = 'translateX(' + offsetX + 'px' + ')' + 
        'translateY(' + offsetY + 'px' + ')' + `scale(${scale})`;

    startPanX = event.clientX;
    startPanY = event.clientY;

    lastScrollPosX = startPanX;
    lastScrollPosY = startPanY;
}

// ========================================================================

function zoom(event) {

    if (switchingToWheelScroll) {
        switchingToWheelScroll = false;
        modelStage.style.transform = 'translateX(' + offsetX + 'px' + ')' +
        'translateY(' + offsetY + 'px' + ')' + `scale(${scale})`;
    }

    let scrollPosX;
    let scrollPosY;

    if (event.target == zoomInButton || event.target == zoomOutButton) {
        scrollPosX = lastScrollPosX;
        scrollPosY = lastScrollPosY;
    } else {
        scrollPosX = event.clientX === undefined ? lastScrollPosX : event.clientX;
        scrollPosY = event.clientY === undefined ? lastScrollPosY : event.clientY;
    }

    if (scrollPosX == 0 && scrollPosY == 0) {
        offsetX = 0;
        offsetY = 0;

        scale += event.deltaY * -0.001;
        scale = Math.min(Math.max(.125, scale), 4);
    } else {
        let stageBoundingBox = modelStage.getBoundingClientRect();
        let posXInStageBeforeZoom = scrollPosX - stageBoundingBox.left;
        let posYInStageBeforeZoom = scrollPosY - stageBoundingBox.top;
        
        let oldScale = scale;
        scale += event.deltaY * -0.001;
        scale = Math.min(Math.max(.125, scale), 4);
    
        let newStageWidthAfterZoom = stageBoundingBox.width * scale/oldScale;
        let newStageHeightAfterZoom = stageBoundingBox.height * scale/oldScale;
    
        let mousePosXPercent = posXInStageBeforeZoom/stageBoundingBox.width;
        let mousePosYPercent = posYInStageBeforeZoom/stageBoundingBox.height;
    
        let posXInStageAfterZoom = mousePosXPercent * newStageWidthAfterZoom;
        let posYInStageAfterZoom = mousePosYPercent * newStageHeightAfterZoom;

        offsetX += posXInStageBeforeZoom - posXInStageAfterZoom;
        offsetY += posYInStageBeforeZoom - posYInStageAfterZoom;
    }

    modelStage.style.transform = 'translateX(' + offsetX + 'px' + ')' +
        'translateY(' + offsetY + 'px' + ')' + `scale(${scale})`;

    lastScrollPosX = event.clientX === undefined ? lastScrollPosX : event.clientX;
    lastScrollPosY = event.clientY === undefined ? lastScrollPosY : event.clientY;

}

// ========================================================================

function getZoomType(delta) {
    return Math.sign(delta) > 0 ? 'in' : 'out'
}

// ========================================================================

function captureMousePosition(event) {
    posX = event.clientX;
    posY = event.clientY;
}

// ========================================================================

function resetView() {
    modelStage.style.transformOrigin = "50% 50% 0";
    lastScrollPosX = 0;
    lastScrollPosY = 0;
    scale = 1;
    posX = 0; posY = 0;
    offsetX = 0; offsetY = 0;
    startPanX = 0; startPanY = 0;
    zeroOriginActive = false;
    modelStage.style.transform = 'translateX(0px) translateY(0px) scale(1)';

    if (window.innerHeight < window.innerWidth) {
        modelImage.style.maxWidth = '95%';
        modelImage.style.height = '95%';
    } else {
        modelImage.style.width = '95%';
        modelImage.style.maxHeight = '95%';
    }
}

// ========================================================================

function wheelZoom(event) {
    if (!zeroOriginActive) {
        let stageBoundingBox = modelStage.getBoundingClientRect();
        offsetX = stageBoundingBox.left;
        offsetY = stageBoundingBox.top;
        modelStage.style.transformOrigin = "0 0";
        zeroOriginActive = true;
        switchingToWheelScroll = true;
    }
    zoom(event);
}

// ========================================================================

function zoomIn(event) {
    zoomInBtnPressed = true;
    zoomInTimer = setInterval(function(){
        if(zoomInBtnPressed) {
            zoom({target: zoomInButton, deltaY: -10});
        } else {
            clearInterval(zoomInTimer);
        }
    }, 25);
}

// ========================================================================

function stopZoomIn() {
    zoomInBtnPressed = false;
}

// ========================================================================

function zoomOut(event) {
    zoomOutBtnPressed = true;
    zoomInTimer = setInterval(function(){
        if(zoomOutBtnPressed) {
            zoom({target: zoomOutButton, deltaY: 10});
        } else {
            clearInterval(zoomInTimer);
        }
    }, 25);
}

// ========================================================================

function stopZoomOut() {
    zoomOutBtnPressed = false;
}

// ========================================================================

function getStageOffset() {

    let style = window.getComputedStyle(modelStage);
    let matrix = new WebKitCSSMatrix(style.transform);

    stageOffsetX = matrix.m41;
    stageOffsetY = matrix.m42;

    return { x: stageOffsetX, y: stageOffsetY };ssss
}

// ========================================================================

document.body.onwheel = wheelZoom;
document.body.addEventListener("mousedown", mousePressed);
document.body.addEventListener("mouseup", disable);
document.body.addEventListener("mouseleave", disable);
document.body.addEventListener("mousemove", pan);

resetButton.addEventListener("click", resetView);
zoomInButton.addEventListener("mousedown", zoomIn);
zoomInButton.addEventListener("mouseup", stopZoomIn);
zoomInButton.addEventListener("mouseleave", stopZoomIn);
zoomOutButton.addEventListener("mousedown", zoomOut);
zoomOutButton.addEventListener("mouseup", stopZoomOut);
zoomOutButton.addEventListener("mouseleave", stopZoomOut);

// ========================================================================
