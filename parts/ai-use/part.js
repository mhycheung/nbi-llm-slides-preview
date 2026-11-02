(function () {
  "use strict";
  // The href is set here because the build refuses any href="http..." in part.html (its
  // check is for files loaded from the network).
  document.querySelectorAll("a.ai-use-link").forEach(function (a) {
    a.href = a.getAttribute("data-url");
  });
})();
