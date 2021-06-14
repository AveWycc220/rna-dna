function render(val) {
    if (typeof val === 'string') {
        document.querySelector('.result-error').innerHTML = val
    } else if (val.length === 3) {
        document.querySelector('.result-error').innerHTML = `${val[0]}`
        document.querySelector('.result-cosine-info').innerHTML = `${val[1]}`
        document.querySelector('.result-cosine-title').innerHTML = `Расстояние городских кварталов`
        document.querySelector('.result-correlation-info').innerHTML = `${val[2]}`
        document.querySelector('.result-correlation-title').innerHTML = `Расстояние Минковского`
    } else {
        document.querySelector('.result-error').innerHTML = ''
        document.querySelector('.result-cosine-title').innerHTML = `Коэффициент Отиаи`
        document.querySelector('.result-correlation-title').innerHTML = `Корреляционный коэффицент`
        document.querySelector('.result-cosine-info').innerHTML = `${val[1]}`
        document.querySelector('.result-correlation-info').innerHTML = `${val[0]}`
    }
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
    await eel.run_programm()().then((data) => render(data)).catch((e) => {
        document.querySelector('.main').innerHTML = e.errorText
    })
}

document.querySelector('.select-1').onclick = selectFileFirst;
document.querySelector('.select-2').onclick = selectFileSecond;
document.querySelector('.run').onclick = run;
