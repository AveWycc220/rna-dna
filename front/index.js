function render(val) {
    console.log(val)
    document.querySelector('.result-info').innerHTML = val
}

function renderPath(val, type) {
    if (type === 1) {
        document.querySelector('.path-1').innerHTML = val
    } else if (type === 2) {
        document.querySelector('.path-2').innerHTML = val
    }
}

async function selectFileFirst() {
    await eel.select_file(1)().then((data) => renderPath(data, 1));
}

async function selectFileSecond() {
    await eel.select_file(2)().then((data) => renderPath(data, 2));
}

async function run() {
    await eel.run_programm()().then((data) => render(data));
}

document.querySelector('.select-1').onclick = selectFileFirst;
document.querySelector('.select-2').onclick = selectFileSecond;
document.querySelector('.run').onclick = run;
