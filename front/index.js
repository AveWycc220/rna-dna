function render(val) {
    if (typeof val === 'string') {
        document.querySelector('.result-error').innerHTML = val
    } else if (val.length === 3) {
        document.querySelector('.result-error').innerHTML = `${val[0]}`
        document.querySelector('.result-cosine-info').innerHTML = `${val[1]}`
        document.querySelector('.result-cosine-title').innerHTML = `Расстояние городских кварталов`
        document.querySelector('.result-correlation-info').innerHTML = `${val[2]}`
        document.querySelector('.result-correlation-title').innerHTML = `Расстояние Минковского`
    } else if (val.length === 2) {
        document.querySelector('.result-error').innerHTML = ''
        document.querySelector('.result-cosine-title').innerHTML = `Коэффициент Отиаи`
        document.querySelector('.result-correlation-title').innerHTML = `Корреляционный коэффицент`
        document.querySelector('.result-cosine-info').innerHTML = `${val[1]}`
        document.querySelector('.result-correlation-info').innerHTML = `${val[0]}`
    }
}

function renderPath(val, type) {
    if (type === 1) {
        document.querySelector('.path-1').innerHTML = val[0]
        document.querySelector('.path-1-length').innerHTML = `Длина: ${val[1]}`
    } else if (type === 2) {
        document.querySelector('.path-2').innerHTML = val[0]
        document.querySelector('.path-2-length').innerHTML = `Длина: ${val[1]}`
    } else if (type === 3) {
        document.querySelector('.path-3').innerHTML = val[0]
    }
}

async function selectFileFirst() {
    await eel.select_file(1)().then((data) => renderPath(data, 1));
}

async function selectFileSecond() {
    await eel.select_file(2)().then((data) => renderPath(data, 2));
}

async function selectFileThird() {
    await eel.select_file(3)().then((data) => renderPath(data, 3));
}

async function showPlot() {
    document.querySelector('.loader-container').classList.remove('disabled')
    document.querySelector('.main').classList.add('disabled')
    document.querySelector('.app').classList.add('loading')
    await eel.show_plot()()
    document.querySelector('.loader-container').classList.add('disabled')
    document.querySelector('.main').classList.remove('disabled')
    document.querySelector('.app').classList.remove('loading')
}

async function showSliding() {
    document.querySelector('.loader-container').classList.remove('disabled')
    document.querySelector('.main').classList.add('disabled')
    document.querySelector('.app').classList.add('loading')
    await eel.open_sliding()()
    document.querySelector('.loader-container').classList.add('disabled')
    document.querySelector('.main').classList.remove('disabled')
    document.querySelector('.app').classList.remove('loading')
}

async function run() {
    document.querySelector('.loader-container').classList.remove('disabled')
    document.querySelector('.main').classList.add('disabled')
    document.querySelector('.app').classList.add('loading')
    if (document.querySelector('.select-control').selectedIndex == 0) {
        await eel.run_programm()().then((data) => render(data)).catch((e) => {
            document.querySelector('.main').innerHTML = e.errorText
        })
    } else {
        const threshold = document.querySelector('.threshold').value
        console.log(threshold)
        if (threshold === '' || (parseFloat(threshold) > 100.0 || parseFloat(threshold) < 0.0)) {
            document.querySelector('.result-error').innerHTML = 'Не введены значения или введены неверно'
        } else {
            await eel.run_programm(parseFloat(threshold))().then((data) => render(data)).catch((e) => {
                document.querySelector('.main').innerHTML = e.errorText
            })
        }
    }
    document.querySelector('.loader-container').classList.add('disabled')
    document.querySelector('.main').classList.remove('disabled')
    document.querySelector('.app').classList.remove('loading')
}

document.querySelector('.select-1').onclick = selectFileFirst;
document.querySelector('.select-2').onclick = selectFileSecond;
document.querySelectorAll('.run')[1].onclick = run;

document.querySelector('.select-control').addEventListener('change', (e) => {
    if (e.target.selectedIndex == 1) {
        document.querySelector('.inputs').classList.remove('disabled')
        document.querySelector('.result').innerHTML = `
                <span class="result-error"></span>
                <button class="open-plot">Открыть и сохранить результаты формулы Бернулли</button>
                <button class="open-sliding">Открыть и сохранить результаты скользящего сравнения</button>
                <button class="open-align">Открыть и сохранить результаты выравнивания</button>
                <button class="open-result">Выбрать папку для записи результатов</button>
                <span class="path-3"></span>
`
        document.querySelector('.open-result').onclick = selectFileThird;
        document.querySelector('.open-plot').onclick = showPlot;
        document.querySelectorAll('.run')[1].classList.add('disabled')
        document.querySelector('.open-sliding').onclick = showSliding;
        document.querySelector('.open-align').onclick = run;
    } else {
        document.querySelector('.inputs').classList.add('disabled')
        document.querySelectorAll('.run')[1].classList.remove('disabled')
        document.querySelector('.result').innerHTML = `
        <span class="result-title">Результат</span>
                <span class="result-error"></span>
                <div class="result-info">
                    <div class="result-cosine">
                        <span class="result-cosine-title">Коэффициент Отиаи</span>
                        <span class="result-cosine-info"></span>
                    </div>
                    <div class="result-correlation">
                        <span class="result-correlation-title">Корреляционный коэффицент</span>
                        <span class="result-correlation-info"></span>
                    </div>
                </div>
`
    }
})
